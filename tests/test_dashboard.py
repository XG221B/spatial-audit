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
        assert page.locator("#experiments tbody tr").count() == 2
        assert page.locator("#experiments a[download]").count() == 6
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


def test_diagnostic_table_sort_search_svg_language_and_mobile(tmp_path):
    import xml.etree.ElementTree as ET

    browser_api = pytest.importorskip("playwright.sync_api")
    with browser_api.sync_playwright() as pw:
        if not Path(pw.chromium.executable_path).is_file():
            pytest.skip("Chromium binary not available")
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 375, "height": 740})
        page.goto((ROOT / "docs/index.html").as_uri())
        page.locator("#kind").select_option("synthetic")
        page.get_by_role("button", name="synthetic-composition", exact=True).click()
        table = page.locator(".feature-diagnostics")
        assert table.locator("tbody tr").count() == 8
        table.get_by_role("button", name="Pooled q").click()
        assert table.locator("th").nth(1).get_attribute("aria-sort") == "ascending"
        sort = table.get_by_role("button", name="Raw energy Δ")
        sort.click()
        ascending = table.locator("tbody td:nth-child(4)").all_text_contents()
        assert list(map(float, ascending)) == sorted(map(float, ascending))
        sort.click()
        assert table.locator("tbody td:nth-child(4)").all_text_contents() == ascending[::-1]
        page.get_by_role("searchbox", name="Search features").fill("F2")
        assert table.locator("tbody tr").count() == 1
        assert table.locator("tbody td").first.inner_text() == "F2"
        with page.expect_download() as pending:
            page.get_by_role("button", name="Export SVG").first.click()
        dest = tmp_path / pending.value.suggested_filename
        pending.value.save_as(dest)
        root = ET.parse(dest).getroot()
        assert root.tag == "{http://www.w3.org/2000/svg}svg" and len(root) > 1
        for expected, heading in [("zh-CN", "逐特征诊断"), ("en", "Per-feature Diagnostics")]:
            page.locator("#language").click()
            assert page.locator("html").get_attribute("lang") == expected
            assert table.get_by_role("heading", name=heading).is_visible()
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        browser.close()
