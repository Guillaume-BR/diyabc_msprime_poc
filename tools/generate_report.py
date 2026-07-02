import subprocess
from pathlib import Path

# Racine du projet = parent de tools/
ROOT = Path(__file__).resolve().parent.parent
NOTES_DIR = ROOT / "notes"
NOTES_DIR.mkdir(exist_ok=True)


def run(cmd, cwd=None):
    return subprocess.check_output(cmd, shell=True, text=True, cwd=cwd or ROOT)


# 1. Tree — lancé depuis la racine
def generate_tree():
    print("Generating tree...")
    out = run("tree -L 3")  # cwd=ROOT, donc l'arbre part bien de la racine
    (NOTES_DIR / "tree.md").write_text(
        f"# Project Tree\n\n```\n{out}\n```", encoding="utf-8"
    )


# 2. Git commits — idem
def generate_commits():
    print("Generating commits...")
    out = run("git log --oneline")  # cwd=ROOT
    (NOTES_DIR / "commits.md").write_text(
        f"# Git Commits\n\n```\n{out}\n```", encoding="utf-8"
    )


# 3. API
def generate_api():
    print("Generating API...")
    run("python3 tools/generate_api_md.py")  # cwd=ROOT
    (ROOT / "API.md").replace(NOTES_DIR / "api.md")


# 4. Report global
def generate_report():
    print("Generating final report...")
    tree = (NOTES_DIR / "tree.md").read_text(encoding="utf-8")
    commits = (NOTES_DIR / "commits.md").read_text(encoding="utf-8")
    api = (NOTES_DIR / "api.md").read_text(encoding="utf-8")

    content = f"""# 📊 Project Report

## 🗂 Structure
{tree}

---

## 📜 Git History
{commits}

---

## 📚 API Documentation
{api}
"""
    (NOTES_DIR / "report.md").write_text(content, encoding="utf-8")


def main():
    generate_tree()
    generate_commits()
    generate_api()
    generate_report()
    print(f"✅ Reports generated in {NOTES_DIR}")


if __name__ == "__main__":
    main()
