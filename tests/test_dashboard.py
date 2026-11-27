"""End-to-end, offline browser checks. Install [browser] and Playwright Chromium."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_file_dashboard_upload_sandbox_and_mobile(tmp_path):
    browser_api = pytest.importorskip("playwright.sync_api")
    with browser_api.sync_playwright() as pw:
        if not Path(pw.chromium.executable_path).is_file():
            pytest.skip("Chromium binary not available")
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1280, "height": 820})
        errors = []
        page.on("pageerror", lambda err: errors.append(str(err)))
        page.goto((ROOT / "docs/index.html").as_uri())
        assert page.locator("#rows tr").count() == 20
        page.locator("#kind").select_option("supplementary")
        assert page.locator("#rows tr").count() == 4
        page.locator("#kind").select_option("variants")
        assert page.locator("#rows tr").count() == 2
        page.locator("#kind").select_option("public")
        page.locator("#search").fill("MTBLS12782")
        assert page.locator("#rows tr").count() == 1
        page.locator("#search").fill("")
        page.locator("#kind").select_option("synthetic")
        assert page.locator("#rows tr").count() == 4
        page.get_by_role("button", name="synthetic-stage", exact=True).click()
        assert page.locator("#detail h2").inner_text() == "synthetic-stage"
        initial = page.locator("#readouts").inner_text()
        page.locator("#n").fill("12")
        assert page.locator("#readouts").inner_text() != initial
        assert page.locator("#curve svg").count() == 3
        raw = ROOT / "results/public/GSE220442.json"
        page.locator("#upload").set_input_files(raw)
        page.wait_for_function(
            "document.querySelector('#notice').textContent.includes('Opened locally')"
        )
        assert page.locator("#kind").input_value() == "local"
        assert page.locator("#rows tr").count() == 1
        assert "Local upload" in page.locator("#detail .tag").inner_text()
        assert page.locator("#detail").get_by_text("Input receipts", exact=False).count() == 0
        page.locator("#kind").select_option("public")
        page.locator("#search").fill("GSE220442")
        page.get_by_role("button", name="GSE220442", exact=True).click()
        assert "Public archive" in page.locator("#detail .tag").inner_text()
        bad = tmp_path / "bad.json"
        bad.write_text('{"unrelated":true}', encoding="utf-8")
        page.locator("#upload").set_input_files(bad)
        page.wait_for_function(
            "document.querySelector('#notice').textContent.includes('Unsupported')"
        )
        page.set_viewport_size({"width": 375, "height": 740})
        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        assert not errors
        browser.close()
