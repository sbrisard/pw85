import argparse
import hashlib
import os
import os.path
import shutil
import urllib.request


def create_parser():
    description = "Download reference data for the PW85 library."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--pw85_ref_data_path",
        required=True,
        help="Where to save the file (string, full path)",
    )
    parser.add_argument("--http_proxy", default="", help="HTTP proxy")
    parser.add_argument("--https_proxy", default="", help="HTTPS proxy")
    return parser


def download_file(url, path, md5, proxies):
    should_download = True
    if os.path.exists(path):
        hash_md5 = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        if hash_md5.hexdigest() == md5:
            should_download = False
    if should_download:
        with urllib.request.urlopen(url) as r, open(path, "wb") as f:
            shutil.copyfileobj(r, f)


if __name__ == "__main__":
    url = "https://zenodo.org/record/3323683/files/pw85_ref_data-20190712.h5"
    md5 = "6009b2a9abde129113612825ee105546"
    parser = create_parser()
    args = parser.parse_args()

    path = args.pw85_ref_data_path
    os.makedirs(os.path.dirname(path), exist_ok=True)

    proxies = {"http": args.http_proxy, "https": args.https_proxy}
    download_file(url, path, md5, proxies)
