import ast
import os
import unittest


class TestNoExternalPrivateCalls(unittest.TestCase):
    def test_no_external_private_calls(self):
        pkg_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "moleditpy",
                "src",
                "moleditpy",
            )
        )
        self.assertTrue(os.path.isdir(pkg_dir), f"Directory not found: {pkg_dir}")

        mismatches = []

        for root, _, files in os.walk(pkg_dir):
            for file in files:
                if not file.endswith(".py"):
                    continue
                path = os.path.join(root, file)
                with open(path, "r", encoding="utf-8") as f:
                    content = f.read()

                try:
                    tree = ast.parse(content, filename=path)
                except SyntaxError as e:
                    mismatches.append(f"{file}: AST parse failed: {e}")
                    continue

                for node in ast.walk(tree):
                    if isinstance(node, ast.Attribute):
                        attr = node.attr
                        # Match single underscore private attributes, excluding dunders
                        if attr.startswith("_") and not attr.startswith("__"):
                            # Check if the base object is 'self'
                            is_self = False
                            if (
                                isinstance(node.value, ast.Name)
                                and node.value.id == "self"
                            ):
                                is_self = True

                            # _cls is an allowed class-level self-reference on the class
                            if not is_self and attr != "_cls":
                                rel_path = os.path.relpath(path, pkg_dir)
                                mismatches.append(
                                    f"{rel_path} line {node.lineno}: private attribute '{attr}' accessed on non-self object"
                                )

        if mismatches:
            self.fail(
                "Found external private calls in main app modules:\n"
                + "\n".join(mismatches)
            )


if __name__ == "__main__":
    unittest.main()
