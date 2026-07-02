import ast
from pathlib import Path

# Racine du projet = parent de tools/
ROOT = Path(__file__).resolve().parent.parent


def format_args(node: ast.FunctionDef):
    args = []
    for a in node.args.args:
        if a.annotation:
            args.append(f"{a.arg}: {ast.unparse(a.annotation)}")
        else:
            args.append(a.arg)
    return ", ".join(args)


def format_return(node: ast.FunctionDef):
    if node.returns:
        return ast.unparse(node.returns)
    return "None"


def extract_functions(file_path: Path):
    tree = ast.parse(file_path.read_text(encoding="utf-8"))
    return [node for node in tree.body if isinstance(node, ast.FunctionDef)]


def generate_markdown(src_dir=None, output=None):
    # Chemins par défaut ancrés sur la racine du projet
    src_dir = Path(src_dir) if src_dir else ROOT / "bridge"
    output = Path(output) if output else ROOT / "API.md"

    out = []
    for file in sorted(src_dir.rglob("*.py")):
        rel = file.relative_to(src_dir)
        functions = extract_functions(file)
        if not functions:
            continue  # on saute les fichiers sans fonctions
        out.append(f"\n# 📄 {rel}\n")
        for fn in functions:
            out.append(f"## `{fn.name}`\n")
            out.append("### Signature\n")
            out.append(
                f"```python\n{fn.name}({format_args(fn)}) -> {format_return(fn)}\n```\n"
            )
            doc = ast.get_docstring(fn)
            if doc:
                out.append("### Description\n")
                out.append(doc.strip() + "\n")
            out.append("---\n")

    output.write_text("\n".join(out), encoding="utf-8")
    print(f"API documentation generated in {output}")


if __name__ == "__main__":
    generate_markdown()
