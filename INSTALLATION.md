# Breizorro Installation Guide

## Basic Installation

```bash
pip install breizorro
# or with uv (faster)
uv pip install breizorro
```

This installs the core masking functionality.

---

## ⚠️ For Catalog Generation: Install Optional Dependencies

**Catalog features require additional packages!**

### Option 1: Catalog Only (Recommended for most users)

```bash
pip install breizorro[catalog]
# or with uv
uv pip install breizorro[catalog]
```

**Includes:**
- `photutils` - For source centroid fitting
- `scikit-image` - For contour detection
- `dask[bag]` - For parallel processing optimization

**Use this for:**
- `--outcatalog` - Source catalog generation
- `--source-fitting` - Different fitting methods (gaussian/centroid/moments/windowed)
- Fast parallel source detection

### Option 2: GUI Only

```bash
pip install breizorro[gui]
```

**Includes:**
- `bokeh` - For interactive visualization

**Use this for:**
- `--gui` - Interactive mask visualization

### Option 3: Everything (Catalog + GUI)

```bash
pip install breizorro[all]
# or with uv
uv pip install breizorro[all]
```

**Includes:**
- All catalog features
- Interactive GUI
- Full feature set

---

## Quick Feature Reference

| Feature | Required Package | Install Command |
|---------|-----------------|-----------------|
| Basic masking | (none) | `pip install breizorro` |
| Source catalog | `[catalog]` | `pip install breizorro[catalog]` |
| Interactive GUI | `[gui]` | `pip install breizorro[gui]` |
| Everything | `[all]` | `pip install breizorro[all]` |

---

## Examples

### Basic Masking (No Optional Deps Required)

```bash
# Simple mask generation
breizorro -r image.fits --outfile mask.fits

# With threshold and dilation
breizorro -r image.fits -t 10 --dilate 3 --outfile mask.fits
```

### Catalog Generation (Requires [catalog])

```bash
# Install first!
pip install breizorro[catalog]

# Generate catalog with default (centroid) fitting
breizorro -r image.fits --outcatalog sources.cat

# Use Gaussian fitting (more accurate, slower)
breizorro -r image.fits --outcatalog sources.cat --source-fitting gaussian

# Use moments fitting (fast alternative)
breizorro -r image.fits --outcatalog sources.cat --source-fitting moments
```

### GUI Visualization (Requires [gui])

```bash
# Install first!
pip install breizorro[gui]

# Open interactive GUI
breizorro -r image.fits --gui
```

---

## Error Messages

### "pip install breizorro[all] to use cataloguing feature"

**Solution:** Install catalog dependencies:
```bash
pip install breizorro[catalog]
# or
pip install breizorro[all]
```

### "Running breizorro gui requires optional dependencies"

**Solution:** Install GUI dependencies:
```bash
pip install breizorro[gui]
# or
pip install breizorro[all]
```

---

## Development Installation

For contributors:

```bash
# Clone repository
git clone https://github.com/ratt-ru/breizorro.git
cd breizorro

# Install with uv (recommended)
uv sync --all-groups

# Or with pip
pip install -e .[all,tests,docs]

# Run tests
pytest breizorro/tests/
```

---

## Verifying Installation

```bash
# Check version
breizorro --help | head -1

# Check if catalog features available
python -c "import photutils, skimage, dask; print('Catalog features: OK')"

# Check if GUI available
python -c "import bokeh; print('GUI features: OK')"
```

---

## Upgrading

```bash
# Upgrade to latest version
pip install --upgrade breizorro[all]

# With uv
uv pip install --upgrade breizorro[all]
```

---

## Uninstallation

```bash
pip uninstall breizorro
```

---

## Summary

**Most users should install:**
```bash
pip install breizorro[catalog]
```

This gives you:
- ✅ All masking features
- ✅ Source catalog generation
- ✅ Multiple fitting methods
- ✅ Dask-optimized parallel processing
- ✅ Fast and reliable

**For everything including GUI:**
```bash
pip install breizorro[all]
```
