# Security Policy

## Supported Versions

| Version | Supported |
| ------- | --------- |
| 3.x (latest) | :white_check_mark: |
| < 3.0 | :x: |

Only the latest release receives security fixes.

## Reporting a Vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

To report a vulnerability, send an email to:

**titech.yoko.hiro@gmail.com**

Include the following in your report:

- A description of the vulnerability and its potential impact
- Steps to reproduce or a proof-of-concept
- The version of MoleditPy affected
- Your operating system and Python version

You can expect an acknowledgement within **5 business days**. If the issue is confirmed, a fix will be released as soon as possible and you will be credited in the release notes (unless you prefer to remain anonymous).

## Scope

MoleditPy is a desktop application that runs locally. The primary attack surface is:

- **File parsing**: Malformed molecular structure files (`.mol`, `.xyz`, `.pme`, etc.) loaded via the UI or drag-and-drop
- **Plugin system**: Plugins execute arbitrary Python code — only install plugins from trusted sources
- **SMILES/SMARTS input**: Strings passed to RDKit for parsing

Out of scope: network attacks, server-side issues (MoleditPy has no network component).
