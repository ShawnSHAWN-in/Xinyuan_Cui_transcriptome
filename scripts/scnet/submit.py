#!/usr/bin/env python3
"""Submit one SCNet container. Success means accepted, not analysis completed."""
import argparse
import base64
import hashlib
import hmac
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, build_opener, HTTPRedirectHandler

AUTH_URL = "https://api.scnet.cn/api/user/v3/tokens"
CENTER_URL = "https://www.scnet.cn/ac/openapi/v2/center"
CREATE_PATH = "/ai/openapi/v2/instance-service/task"


class ScnetError(RuntimeError):
    pass


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def required(env, name):
    value = env.get(name, "").strip()
    if not value:
        raise ScnetError(f"Missing configuration: {name}")
    return value


def signature(user, access_key, secret_key, timestamp):
    message = json.dumps(
        {"accessKey": access_key, "timestamp": timestamp, "user": user},
        ensure_ascii=False, separators=(",", ":"), sort_keys=True,
    )
    return hmac.new(secret_key.encode(), message.encode(), hashlib.sha256).hexdigest()


def api(method, url, headers, payload=None):
    # No automatic retry of mutations: a timeout may follow successful creation.
    if urlsplit(url).scheme != "https":
        raise ScnetError("API URL must use HTTPS")
    data = None if payload is None else json.dumps(payload, ensure_ascii=False).encode()
    request = Request(url, data=data, method=method, headers={"Content-Type": "application/json", **headers})
    try:
        with build_opener(NoRedirect).open(request, timeout=60) as response:
            result = json.load(response)
    except HTTPError as exc:
        raise ScnetError(f"API HTTP {exc.code}; response body omitted") from None
    except (URLError, TimeoutError, OSError, ValueError):
        raise ScnetError("API transport/JSON error; verify remote state before retrying") from None
    if not isinstance(result, dict) or str(result.get("code")) != "0":
        code = result.get("code", "missing") if isinstance(result, dict) else "invalid-response"
        raise ScnetError(f"API rejected request, code={code}; response body omitted")
    return result.get("data")


def mask(value):
    # Mask every regional token, including tokens that will not be used.
    print(f"::add-mask::{value}")


def container_url(base_url, path=CREATE_PATH):
    # Center examples include /ai; do not generate /ai/ai/openapi/...
    base_url = base_url.rstrip("/")
    return base_url + (path[3:] if base_url.endswith("/ai") else path)


def login(env):
    user = required(env, "SCNET_USER")
    access_key = required(env, "SCNET_ACCESS_KEY")
    secret_key = required(env, "SCNET_SECRET_KEY")
    timestamp = str(int(time.time()))
    digest = signature(user, access_key, secret_key, timestamp)
    mask(digest)
    regions = api("POST", AUTH_URL, {
        "user": user, "accessKey": access_key, "timestamp": timestamp, "signature": digest,
    }, {})
    if not isinstance(regions, list):
        raise ScnetError("Unexpected authentication response")
    for region in regions:
        if region.get("token"):
            mask(region["token"])
    cluster_id = required(env, "SCNET_CLUSTER_ID")
    matches = [r for r in regions if str(r.get("clusterId")) == cluster_id and r.get("token")]
    if cluster_id == "0" or len(matches) != 1:
        raise ScnetError("Select one enabled compute region; clusterId=0 is a platform token")
    token = matches[0]["token"]
    center = api("GET", CENTER_URL, {"token": token})
    if isinstance(center, list):
        centers = [c for c in center if str(c.get("id")) == cluster_id]
        if len(centers) != 1:
            raise ScnetError("Region lookup was ambiguous")
        center = centers[0]
    if not isinstance(center, dict) or str(center.get("id")) != cluster_id:
        raise ScnetError("Region lookup returned a different region")
    urls = [u["url"] for u in center.get("aiUrls", [])
            if str(u.get("enable")).lower() == "true" and u.get("url")]
    if not urls:
        raise ScnetError("Selected region has no enabled container API")
    base_url = urls[0].rstrip("/")
    if urlsplit(base_url).scheme != "https":
        raise ScnetError("Regional container API must use HTTPS")
    return token, base_url


def build_payload(env):
    config = json.loads(required(env, "SCNET_CONTAINER_JSON"))
    if not isinstance(config, dict) or "REPLACE_WITH" in json.dumps(config):
        raise ScnetError("Replace example resource/image/mount values")
    for key in ("taskType", "acceleratorType", "version", "imagePath", "resourceGroup", "timeoutLimit"):
        if not isinstance(config.get(key), str) or not config[key].strip():
            raise ScnetError(f"Missing container field: {key}")
    for key in ("cpuNumber", "ramSize", "gpuNumber", "taskNumber"):
        if type(config.get(key)) is not int:
            raise ScnetError(f"Container field must be integer: {key}")
    if config["cpuNumber"] < 1 or config["ramSize"] < 1 or config["gpuNumber"] < 0 or config["taskNumber"] != 1:
        raise ScnetError("Require positive CPU/RAM, nonnegative GPU and taskNumber=1")
    if config["acceleratorType"] not in ("cpu", "gpu", "dcu", "mlu") or config["taskType"] not in ("ssh", "jupyter", "codeserver", "rstudio"):
        raise ScnetError("Invalid accelerator/task type")
    limit = re.fullmatch(r"(\d{2}):([0-5]\d):([0-5]\d)", config["timeoutLimit"])
    if not limit or not any(int(n) for n in limit.groups()):
        raise ScnetError("Use a finite nonzero timeoutLimit in HH:MM:SS")
    root = required(env, "SCNET_WORK_ROOT")
    if not root.startswith("/") or ".." in PurePosixPath(root).parts or root == "/":
        raise ScnetError("SCNET_WORK_ROOT must be an absolute mounted directory")
    mounts = config.get("mountInfoList", [])
    if not mounts or not any(
        isinstance(m.get("targetPath"), str) and m["targetPath"].startswith("/")
        and (PurePosixPath(root) == PurePosixPath(m["targetPath"])
             or PurePosixPath(m["targetPath"]) in PurePosixPath(root).parents)
        for m in mounts
    ):
        raise ScnetError("SCNET_WORK_ROOT must be inside a persistent mount")
    for mount in mounts:
        if not mount.get("sourcePath", "").startswith("/") or mount.get("type") not in ("data", "path"):
            raise ScnetError("Invalid mount sourcePath/type")
    repo = required(env, "GITHUB_REPOSITORY")
    sha = required(env, "GITHUB_SHA")
    run_key = required(env, "GITHUB_RUN_ID") + "-" + required(env, "GITHUB_RUN_ATTEMPT")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo) or not re.fullmatch(r"[a-f0-9]{40}", sha):
        raise ScnetError("Invalid repository/commit")
    if not re.fullmatch(r"\d+-\d+", run_key):
        raise ScnetError("Invalid run ID/attempt")
    variables = {
        "SCNET_WORK_ROOT": root,
        "SCNET_RUN_KEY": run_key,
        "SCNET_COMMIT": sha,
        "SCNET_REPO_URL": f"https://github.com/{repo}.git",
        "SCNET_RUN_COMMAND": required(env, "SCNET_RUN_COMMAND"),
        "GIT_TERMINAL_PROMPT": "0",
    }
    entry = Path(__file__).with_name("container_entry.sh").read_text()
    exports = "\n".join(f"export {k}={shlex.quote(v)}" for k, v in variables.items())
    encoded = base64.b64encode((exports + "\n" + entry).encode()).decode()
    config.update({
        "instanceServiceName": f"gh-{run_key}",
        "description": f"{repo}@{sha}; GitHub run {run_key}",
        "useStartScript": True,
        "startScriptActionScope": "header",
        "startScriptContent": f"printf '%s' '{encoded}' | base64 -d | bash\n",
    })
    return config, run_key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="Validate without authentication or API calls")
    args = parser.parse_args()
    payload, run_key = build_payload(os.environ)
    if args.dry_run:
        print(f"Configuration valid; run={run_key}; no API calls made")
        return
    token, base_url = login(os.environ)
    # A lost response may still have created a container. Do not blindly rerun.
    print(f"Submitting instance name {payload['instanceServiceName']}", flush=True)
    instance_id = api("POST", container_url(base_url), {"token": token}, payload)
    if not isinstance(instance_id, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", instance_id):
        raise ScnetError("Unexpected instance ID; verify remote state before retrying")
    print(f"Accepted instance: {instance_id}")
    message = (
        f"SCNet accepted container `{instance_id}` (name `{payload['instanceServiceName']}`).\n\n"
        f"Commit: `{os.environ['GITHUB_SHA']}`\n\n"
        f"Remote log: `{os.environ['SCNET_WORK_ROOT']}/runs/{run_key}/run.log`\n\n"
        "This workflow confirms submission only. Check remote `exit_code`: 0 = success, "
        "nonzero = failed, absent = queued/running/interrupted.\n\n"
        f"Container auto-stop limit: `{payload['timeoutLimit']}`. Script completion does not stop the instance early.\n"
    )
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as handle:
            handle.write(message)


if __name__ == "__main__":
    try:
        main()
    except ScnetError as exc:
        print(f"SCNet submission failed: {exc}", file=sys.stderr)
        sys.exit(1)
    except (ValueError, KeyError, TypeError):
        # Do not echo configuration, tokens, request payloads or response bodies.
        print("SCNet submission failed. Check configuration/API access; if a create request was sent, "
              "check the instance name in SCNet before retrying.", file=sys.stderr)
        sys.exit(1)
