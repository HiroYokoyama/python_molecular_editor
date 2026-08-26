# Project Philosophy

This document describes the design philosophy of MoleditPy.
It is intended for users and contributors who want to understand
the long-term direction and architectural decisions of the project.

## 1. Vision

**"Simplicity at the Core. Limitless Potential in the Extensions."**

We believe a molecular editor should be easy to pick up, yet powerful enough to handle complex research. MoleditPy provides a streamlined foundation for visualization that puts simplicity first, while offering an open architecture that allows for infinite customization and expansion.

## 2. Core Values

### Simplicity First
We believe in keeping the core software simple and fast. Instead of packing every possible feature into the main program, we provide a lightweight foundation. This ensures a responsive experience and allows you to add exactly the features you need via plugins.

### Open and Transparent
Trust is essential for science. MoleditPy is designed to be transparent, allowing you to see exactly how your data is processed. Whether you are checking a calculation or learning how the code works, the logic is always accessible.

### Private by Default
Research data is often confidential long before it is published. The core application never connects to the network: it contains no networking code, sends no telemetry, and checks nothing online. Everything you sketch, convert, and export stays on your machine. Connectivity enters only through plugins you choose to install yourself — the Plugin Installer, for instance, downloads from the plugin registry — so going online is always something you opt into. The one network-facing element in the core is the "Explore Plugins Online" button, which hands a URL to your own web browser and transmits none of your data.

### Freedom to Create
Your research is unique, and your tools should be too. MoleditPy gives you the freedom to go beyond standard menus. Start with intuitive built-in tools, and when you are ready, extend your workflow with Python scripts to automate tasks and interact with molecules in ways that fit your specific goals.

### Connected Ecosystem
MoleditPy works harmoniously with the tools you already use. It seamlessly connects with the Python scientific ecosystem (like RDKit and PyVista), creating a unified environment where visualization and analysis happen together.
