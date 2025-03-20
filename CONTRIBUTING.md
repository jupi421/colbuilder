# Contributing to ColBuilder

Thank you for your interest in contributing to ColBuilder! This guide outlines the basics of how to contribute.

## Quick Start

1. **Fork and clone the repository**
   ```bash
   git clone git@github.com:YOUR-USERNAME/colbuilder.git
   cd colbuilder
   ```

2. **Set up your environment**
   ```bash
   conda create -n colbuilder-dev python=3.9
   conda activate colbuilder-dev
   pip install -e .
   ```

3. **Create a branch and make changes**
   ```bash
   git checkout -b feature/your-feature-name
   # Make your changes
   git commit -m "Add feature: description"
   git push origin feature/your-feature-name
   ```

4. **Open a pull request** from your fork to the main repository

## Reporting Issues

When reporting bugs or requesting features:
- Search existing issues first
- Use descriptive titles
- Include detailed descriptions
- For bugs: add steps to reproduce, expected vs. actual behavior, and system info

## Development Guidelines

- **Code Organization**: 
  - Follow the existing project structure
  - Place new modules in the appropriate directories
  - Maintain separation between core functions and utilities

- **Documentation**: 
  - Add docstrings to all new functions and classes
  - Update README and other documentation as needed
  - Include example usage for new features

## Pull Request Guidelines

- Focus each PR on a single feature or fix
- Provide a clear description of the changes
- Reference related issues
- Be responsive to feedback

## Versioning

ColBuilder follows a simple versioning system:
- Version numbers are in the format X.Y
- Major changes increase X
- Minor improvements increase Y

---

Thank you for helping improve ColBuilder for the scientific community!