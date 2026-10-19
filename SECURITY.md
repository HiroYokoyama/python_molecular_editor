# Security Policy

## Supported Versions

| Version | Supported |
| ------- | --------- |
| 5.x (latest stable release) | :white_check_mark: |
| 4.x and older | :x: |

Only the latest release receives security fixes.

## Reporting a Vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

To report a vulnerability, send an email to:

**titech.yoko.hiro[at]gmail.com**

Include the following in your report:

- A description of the vulnerability and its potential impact
- Steps to reproduce or a proof-of-concept
- The version of MoleditPy affected
- Your operating system and Python version

You can expect an acknowledgement within **5 business days**. If the issue is confirmed, a fix will be released as soon as possible and you will be credited in the release notes (unless you prefer to remain anonymous).

## Scope

MoleditPy is a desktop application that runs locally. The primary attack surface is:

- **File parsing**: Malformed molecular structure files (`.mol`, `.xyz`, `.pmeprj`, `.pmeraw` etc.) loaded via the UI or drag-and-drop
- **Plugin system**: Plugins execute arbitrary Python code — only install plugins from trusted sources
- **SMILES/SMARTS input**: Strings passed to RDKit for parsing
- **Native libraries**: PySide6/Qt, RDKit and VTK/PyVista process data in the application process; include their versions when reporting crashes
- **Plugin installation and integrations**: Downloaded Python code, external programs and any HTTP/MCP/SSH services enabled by installed plugins

The core editor does not start an inbound network service by default. Installed plugins can access files, execute programs and communicate over the network with the user's privileges; they are not sandboxed. Report collection or standalone-plugin issues to the relevant repository, and include the plugin name/version and enabled integrations.

A registry SHA-256 check verifies that a download matches the catalogued bytes; it does not establish that the code is safe. Review the source and publisher before installation. Historical version tags are provided for reproducibility and do not receive security fixes. Third-party runtime vulnerabilities may require updating the affected dependency as well as MoleditPy.
