import sys

if __name__ == '__main__':
    version_mask = '@project_version@'
    with open(sys.argv[1]) as fi:
        version = next(fi).strip()
    with open(sys.argv[2], 'r') as fi:
        lines = [line.replace(version_mask, version) for line in fi]
    with open(sys.argv[3], 'w') as fo:
        fo.writelines(lines)
