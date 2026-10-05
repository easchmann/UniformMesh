import configparser
import sys


def define(name, section):
    head = name
    if 'args' in section:
        args = [a.strip() for a in section['args'].split(',')]
        head = f"{name}({','.join(args)})"

    lines = section['definition'].strip().splitlines()
    if len(lines) <= 1:
        return f"#define {head} {lines[0] if lines else ''}"
    body = " \\\n    ".join(line.rstrip() for line in lines)
    return f"#define {head} {body}"


def main():
    args = sys.argv[1:]
    out_path = 'mesh_config.h'
    if '-o' in args:
        i = args.index('-o')
        out_path = args[i + 1]
        args = args[:i] + args[i + 2:]

    config = configparser.ConfigParser(
        interpolation=None, delimiters=('=',), comment_prefixes=(';',))
    config.optionxform = str
    for path in args:
        config.read(path)

    with open(out_path, 'w') as out:
        out.write("#ifndef MESH_CONFIG_H\n#define MESH_CONFIG_H\n\n")
        for name in config.sections():
            out.write(define(name, config[name]) + "\n\n")
        out.write("#endif\n")


if __name__ == '__main__':
    main()
