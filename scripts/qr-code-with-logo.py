"""Generate OPAL-branded QR codes with a white-backed center logo."""

from pathlib import Path
import sys

import qrcode
from PIL import Image


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LOGO_PATH = REPOSITORY_ROOT / "assets" / "static" / "opal-HEX-transparent.png"
DEFAULT_LOGO_SCALE = 0.18


def add_center_logo(qr_image, logo_path=DEFAULT_LOGO_PATH, logo_scale=DEFAULT_LOGO_SCALE):
    """Overlay a centered transparent logo on a padded, white backing."""
    if not 0 < logo_scale < 1:
        raise ValueError("logo_scale must be greater than 0 and less than 1.")

    logo_path = Path(logo_path)
    if not logo_path.is_file():
        raise FileNotFoundError(f"Logo file not found: {logo_path}")

    qr_image = qr_image.convert("RGBA")
    logo = Image.open(logo_path).convert("RGBA")

    qr_width, qr_height = qr_image.size
    max_logo_dimension = int(min(qr_width, qr_height) * logo_scale)
    logo.thumbnail((max_logo_dimension, max_logo_dimension), Image.Resampling.LANCZOS)

    padding = max(4, max_logo_dimension // 10)
    backing = Image.new(
        "RGBA",
        (logo.width + 2 * padding, logo.height + 2 * padding),
        "white",
    )
    logo_position = (
        (backing.width - logo.width) // 2,
        (backing.height - logo.height) // 2,
    )
    backing.alpha_composite(logo, logo_position)

    backing_position = (
        (qr_width - backing.width) // 2,
        (qr_height - backing.height) // 2,
    )
    qr_image.alpha_composite(backing, backing_position)

    return qr_image


def generate_qr_code(
    url,
    filename=None,
    size=10,
    border=4,
    logo_path=DEFAULT_LOGO_PATH,
    logo_scale=DEFAULT_LOGO_SCALE,
):
    """Generate an OPAL-branded QR code for a URL.

    Args:
        url (str): The website URL to encode.
        filename (str | Path): Optional output PNG filename.
        size (int): Size of each QR box in pixels.
        border (int): Quiet-zone border size in boxes.
        logo_path (str | Path): Transparent logo placed over the QR center.
        logo_scale (float): Maximum logo dimension relative to QR width.
    """
    if not url.startswith(("http://", "https://")):
        url = "https://" + url

    qr = qrcode.QRCode(
        version=None,
        error_correction=qrcode.constants.ERROR_CORRECT_H,
        box_size=size,
        border=border,
    )
    qr.add_data(url)
    qr.make(fit=True)

    image = qr.make_image(fill_color="black", back_color="white")
    image = add_center_logo(image, logo_path, logo_scale)

    if filename is None:
        safe_url = (
            url.replace("https://", "")
            .replace("http://", "")
            .replace("/", "_")
            .replace(":", "")
        )
        filename = f"qr_{safe_url}.png"

    output_path = Path(filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(output_path)
    print(f"QR code saved as: {output_path}")
    print(f"QR code contains: {url}")

    return output_path


def batch_generate_qr_codes(urls, output_dir="qr_codes"):
    """Generate white-backed OPAL QR codes for multiple URLs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_files = []
    for index, url in enumerate(urls, start=1):
        safe_url = (
            url.replace("https://", "")
            .replace("http://", "")
            .replace("/", "_")
            .replace(":", "")
        )
        filename = output_dir / f"qr_{index}_{safe_url}.png"

        try:
            generated_files.append(generate_qr_code(url, filename))
        except Exception as error:
            print(f"Error generating QR code for {url}: {error}")

    return generated_files


def main():
    """Provide an interactive command-line interface."""
    print("OPAL QR Code Generator")
    print("=" * 23)

    while True:
        print("\nOptions:")
        print("1. Generate single QR code")
        print("2. Generate multiple QR codes")
        print("3. Exit")

        choice = input("\nEnter your choice (1-3): ").strip()

        if choice == "1":
            url = input("Enter website URL: ").strip()
            if not url:
                print("Please enter a valid URL.")
                continue

            try:
                filename = input("Enter filename (press Enter for default): ").strip()
                generate_qr_code(url, filename or None)
            except Exception as error:
                print(f"Error: {error}")

        elif choice == "2":
            print("Enter URLs (one per line, empty line to finish):")
            urls = []
            while url := input().strip():
                urls.append(url)

            if not urls:
                print("No URLs provided.")
                continue

            try:
                output_dir = input("Enter output directory (press Enter for 'qr_codes'): ").strip()
                generated_files = batch_generate_qr_codes(urls, output_dir or "qr_codes")
                print(f"\nGenerated {len(generated_files)} QR codes in '{output_dir or 'qr_codes'}' directory")
            except Exception as error:
                print(f"Error: {error}")

        elif choice == "3":
            print("Goodbye!")
            break

        else:
            print("Invalid choice. Please enter 1, 2, or 3.")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        url = sys.argv[1]
        filename = sys.argv[2] if len(sys.argv) > 2 else None
        generate_qr_code(url, filename)
    else:
        main()

# Example usage:
# py scripts/qr-code-with-logo.py https://connect.fisheries.noaa.gov/opal/ assets/images/test-qr-logo-opal-home.png
