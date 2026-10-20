import base64
from contextlib import redirect_stdout
import hashlib
import hmac
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("scnet_submit", ROOT / "scripts/scnet/submit.py")
scnet = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(scnet)


def environment(root="/mnt/work"):
    config = json.loads((ROOT / "scripts/scnet/container.example.json").read_text())
    config.update(version="r:4", imagePath="registry/r:4", resourceGroup="cpu")
    config["mountInfoList"][0]["sourcePath"] = "/storage/user/project"
    return {
        "SCNET_CONTAINER_JSON": json.dumps(config), "SCNET_WORK_ROOT": root,
        "SCNET_RUN_COMMAND": "Rscript -e 'targets::tar_make()'",
        "GITHUB_REPOSITORY": "owner/repo", "GITHUB_SHA": "a" * 40,
        "GITHUB_RUN_ID": "123", "GITHUB_RUN_ATTEMPT": "1",
        "SCNET_CLUSTER_ID": "12", "SCNET_USER": "user",
        "SCNET_ACCESS_KEY": "ak", "SCNET_SECRET_KEY": "sk",
    }


class SubmissionTests(unittest.TestCase):
    def test_signature_official_canonical_order(self):
        expected = hmac.new(b"secret", b'{"accessKey":"ak","timestamp":"100","user":"bob"}', hashlib.sha256).hexdigest()
        self.assertEqual(scnet.signature("bob", "ak", "secret", "100"), expected)

    def test_api_prefix(self):
        expected = "https://center.example/ai/openapi/v2/instance-service/task"
        self.assertEqual(scnet.container_url("https://center.example/ai/"), expected)
        self.assertEqual(scnet.container_url("https://center.example"), expected)

    def test_enabled_regional_token(self):
        responses = [
            [{"clusterId": "0", "token": "platform"}, {"clusterId": "12", "token": "regional"}],
            {"id": 12, "aiUrls": [{"enable": False, "url": "https://disabled"},
                                   {"enable": "true", "url": "https://center/ai"}]},
        ]
        with patch.object(scnet, "api", side_effect=responses) as api, redirect_stdout(io.StringIO()):
            token, url = scnet.login(environment())
        self.assertEqual((token, url), ("regional", "https://center/ai"))
        self.assertEqual(api.call_args_list[1].args[2], {"token": "regional"})

    def test_payload_and_shell_quoting(self):
        env = environment("/mnt/work/带 空格'目录")
        env["SCNET_RUN_COMMAND"] = "printf '%s' '$(touch /should-not-run)'"
        payload, key = scnet.build_payload(env)
        self.assertEqual(key, "123-1")
        self.assertTrue(payload["useStartScript"])
        self.assertEqual(payload["startScriptActionScope"], "header")
        encoded = payload["startScriptContent"].split("'")[3]
        script = base64.b64decode(encoded).decode()
        subprocess.run(["bash", "-n"], input=script, text=True, check=True)
        exports = script.split("#!/usr/bin/env bash")[0]
        result = subprocess.run(["bash", "-c", exports + "\nprintf '%s' \"$SCNET_RUN_COMMAND\""], text=True, capture_output=True, check=True)
        self.assertEqual(result.stdout, env["SCNET_RUN_COMMAND"])

    def test_reject_unbounded_and_unmounted(self):
        env = environment()
        config = json.loads(env["SCNET_CONTAINER_JSON"])
        config["timeoutLimit"] = "unlimited"
        env["SCNET_CONTAINER_JSON"] = json.dumps(config)
        with self.assertRaises(scnet.ScnetError):
            scnet.build_payload(env)
        with self.assertRaises(scnet.ScnetError):
            scnet.build_payload(environment("/outside/mount"))

    def test_dry_run_has_no_api_calls(self):
        with patch.dict(os.environ, environment(), clear=True), patch.object(scnet, "api") as api:
            with patch("sys.argv", ["submit.py", "--dry-run"]), redirect_stdout(io.StringIO()):
                scnet.main()
        api.assert_not_called()

    def test_http_success_api_failure(self):
        response = unittest.mock.MagicMock()
        response.__enter__.return_value.read.return_value = b'{"code":"10008","msg":"sensitive"}'
        with patch.object(scnet, "build_opener") as opener:
            opener.return_value.open.return_value = response
            with self.assertRaisesRegex(scnet.ScnetError, "code=10008"):
                scnet.api("POST", "https://center/api", {}, {})


class RemoteEntryTests(unittest.TestCase):
    def test_exact_commit_cache_reuse_and_failure_status(self):
        with tempfile.TemporaryDirectory(prefix="scnet test ") as temp:
            folder = Path(temp)
            origin = folder / "origin"
            origin.mkdir()

            def git(*args):
                return subprocess.run(["git", "-C", str(origin), *args], check=True, text=True, capture_output=True).stdout.strip()

            git("init")
            (origin / "version.txt").write_text("first")
            git("add", ".")
            git("-c", "user.name=Test", "-c", "user.email=test@example.org", "commit", "-m", "first")
            first = git("rev-parse", "HEAD")
            (origin / "version.txt").write_text("second")
            git("add", ".")
            git("-c", "user.name=Test", "-c", "user.email=test@example.org", "commit", "-m", "second")
            work = folder / "persistent mount"
            env = {**os.environ, "SCNET_WORK_ROOT": str(work), "SCNET_RUN_KEY": "123-1",
                   "SCNET_COMMIT": first, "SCNET_REPO_URL": str(origin),
                   "SCNET_RUN_COMMAND": "mkdir -p _targets; cp version.txt _targets/cached"}
            script = ROOT / "scripts/scnet/container_entry.sh"
            subprocess.run(["bash", str(script)], env=env, check=True)
            self.assertEqual((work / "project/version.txt").read_text(), "first")
            self.assertEqual((work / "runs/123-1/exit_code").read_text().strip(), "0")
            env.update(SCNET_RUN_KEY="123-2", SCNET_RUN_COMMAND="test -f _targets/cached; exit 17")
            result = subprocess.run(["bash", str(script)], env=env)
            self.assertEqual(result.returncode, 17)
            self.assertEqual((work / "runs/123-2/exit_code").read_text().strip(), "17")
            self.assertEqual((work / "project/_targets/cached").read_text(), "first")


if __name__ == "__main__":
    unittest.main()
