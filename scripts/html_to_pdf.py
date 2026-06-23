"""Render a MathJax-enabled Jupyter HTML export to an A4 PDF.

Usage:
    D:\\Anaconda3\\python.exe html_to_pdf.py "Solutions 771-780.html"
"""

from __future__ import annotations

import argparse
import base64
import os
import sys
from pathlib import Path

from selenium import webdriver
from selenium.common.exceptions import TimeoutException, WebDriverException
from selenium.webdriver.common.print_page_options import PrintOptions
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.ui import WebDriverWait

DEFAULT_FIREFOX = Path(r"C:\Program Files\Mozilla Firefox\firefox.exe")
DEFAULT_GECKODRIVER = Path(r"C:\Program Files\Mozilla Firefox\geckodriver.exe")

PRINT_STYLE = r"""
@media print {
    .container, div#notebook {
        width: 100% !important;
        max-width: none !important;
        margin: 0 !important;
        padding: 0 !important;
    }
    body { margin: 0 !important; padding: 0 !important; }

    /* A shared prompt column keeps Markdown aligned with code and output:
       In [12]: <code>
                <Markdown> */
    div.input, div.output_area, div.text_cell {
        width: 100% !important;
        display: flex !important;
        flex-direction: row !important;
        align-items: stretch !important;
    }
    div.prompt, div.input_prompt, div.output_prompt {
        flex: 0 0 50px !important;
        min-width: 50px !important;
        font-size: 9pt !important;
        line-height: 1.1 !important;
        text-align: left !important;
        padding: 0 0.4em 0 0 !important;
        box-sizing: border-box !important;
    }
    div.inner_cell, div.output_subarea {
        flex: 1 1 0 !important;
        width: auto !important;
        min-width: 0 !important;
    }
    div.output_subarea {
        max-width: calc(100% - 60px) !important;
    }

    /* Do not print Bootstrap's appended '(URL)' text after links. */
    a[href]::after { content: none !important; }

    div.text_cell_render { font-size: 9pt !important; line-height: 1.0 !important; }
    div.input_area pre, div.output_area pre, .highlight pre {
        font-size: 9pt !important;
        line-height: 1.0 !important;
        white-space: pre-wrap !important;
        overflow-wrap: anywhere !important;
    }

    /* Override Jupyter's 'avoid page break' print rules. */
    div.cell, div.code_cell, div.input, div.output_wrapper, div.output,
    div.output_area, div.input_area, .highlight, pre, .MathJax_Display,
    .MathJax_Display > span, .MathJax_Display .MathJax {
        page-break-inside: auto !important;
        break-inside: auto !important;
    }
}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render a MathJax-enabled Jupyter HTML export as an A4 PDF."
    )
    parser.add_argument("html", type=Path, help="Input HTML file")
    parser.add_argument(
        "--firefox",
        type=Path,
        default=Path(os.environ.get("FIREFOX_BINARY", DEFAULT_FIREFOX)),
        help="Path to firefox.exe (or set FIREFOX_BINARY)",
    )
    parser.add_argument(
        "--geckodriver",
        type=Path,
        default=Path(os.environ.get("GECKODRIVER", DEFAULT_GECKODRIVER)),
        help="Path to geckodriver.exe (or set GECKODRIVER)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=90,
        help="Seconds to wait for MathJax (default: 90)",
    )
    return parser.parse_args()


def wait_for_mathjax(driver: webdriver.Firefox, timeout: int) -> None:
    wait = WebDriverWait(driver, timeout)
    wait.until(lambda d: d.execute_script("return document.readyState === 'complete'"))
    wait.until(
        lambda d: d.execute_script(
            "return Boolean(window.MathJax && window.MathJax.Hub)"
        )
    )
    driver.set_script_timeout(timeout)
    driver.execute_async_script("""
        const done = arguments[0];
        MathJax.Hub.Queue(["Typeset", MathJax.Hub], () => done(true));
    """)


def convert_literal_markdown_links(driver: webdriver.Firefox) -> None:
    """Convert prose [label](https://example.com) text to anchors."""
    driver.execute_script(r"""
        const pattern = new RegExp(
            "\\[([^\\]\\r\\n]+)\\]\\((https?:\\/\\/[^\\s)]+)\\)",
            "g"
        );
        const root = document.querySelector("#notebook") || document.body;
        const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
            acceptNode(node) {
                pattern.lastIndex = 0;
                if (!pattern.test(node.nodeValue)) return NodeFilter.FILTER_REJECT;

                const parent = node.parentElement;
                if (!parent || parent.closest(
                    "a, pre, code, script, style, .MathJax, .MathJax_Display"
                )) return NodeFilter.FILTER_REJECT;

                return NodeFilter.FILTER_ACCEPT;
            }
        });

        const nodes = [];
        while (walker.nextNode()) nodes.push(walker.currentNode);

        for (const node of nodes) {
            const text = node.nodeValue;
            const fragment = document.createDocumentFragment();
            let last = 0;
            let match;
            pattern.lastIndex = 0;

            while ((match = pattern.exec(text)) !== null) {
                fragment.append(text.slice(last, match.index));
                const link = document.createElement("a");
                link.href = match[2];
                link.textContent = match[1];
                fragment.append(link);
                last = pattern.lastIndex;
            }
            fragment.append(text.slice(last));
            node.replaceWith(fragment);
        }
    """)


def inject_print_style(driver: webdriver.Firefox) -> None:
    driver.execute_script(
        """
        const style = document.createElement("style");
        style.textContent = arguments[0];
        document.head.appendChild(style);
        """,
        PRINT_STYLE,
    )


def make_print_options() -> PrintOptions:
    options = PrintOptions()
    options.set_page_size(PrintOptions.A4)
    options.margin_top = 1.3
    options.margin_bottom = 1.3
    options.margin_left = 1.3
    options.margin_right = 1.3
    options.scale = 1.0
    options.background = True
    return options


def main() -> int:
    args = parse_args()
    html = args.html.resolve()

    if not html.is_file():
        print(f"Input HTML file not found: {html}", file=sys.stderr)
        return 2
    if not args.firefox.is_file():
        print(f"Firefox executable not found: {args.firefox}", file=sys.stderr)
        return 2
    if not args.geckodriver.is_file():
        print(f"Geckodriver executable not found: {args.geckodriver}", file=sys.stderr)
        return 2

    output = html.with_suffix(".pdf")
    firefox_options = Options()
    firefox_options.binary_location = str(args.firefox)
    firefox_options.add_argument("-headless")
    driver: webdriver.Firefox | None = None

    try:
        driver = webdriver.Firefox(
            options=firefox_options,
            service=Service(executable_path=str(args.geckodriver)),
        )
        driver.get(html.as_uri())
        wait_for_mathjax(driver, args.timeout)
        convert_literal_markdown_links(driver)
        inject_print_style(driver)
        output.write_bytes(base64.b64decode(driver.print_page(make_print_options())))
    except TimeoutException:
        print(
            "Timed out while waiting for MathJax. Check the HTML's MathJax URL.",
            file=sys.stderr,
        )
        return 1
    except WebDriverException as error:
        print(f"Firefox PDF rendering failed: {error.msg}", file=sys.stderr)
        return 1
    finally:
        if driver is not None:
            driver.quit()

    print(f"Created: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
