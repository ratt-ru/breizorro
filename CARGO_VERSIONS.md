# Using Different Breizorro Versions with Stimela 2

This document explains how to use different versions of breizorro with Stimela 2.


## Available Versions

### Stable Releases

- **`latest`** (default): Most recent stable release with all features
  ```bash
  stimela run -S breizorro::breizorro.yaml [params...]
  ```

- **`0.2.1`**: Specific stable version
  ```yaml
  # version-override.yaml
  cabs:
    breizorro:
      image:
        version: '0.2.1'
  ```
  
  ```bash
  stimela run -S breizorro::breizorro.yaml version-override.yaml [params...]
  ```

### Development Versions

- **`updates`**: Latest from the updates branch
  ```yaml
  cabs:
    breizorro:
      image:
        version: updates
  ```

- **`main`**: Latest from the main branch
  ```yaml
  cabs:
    breizorro:
      image:
        version: main
  ```

## Override Image Version

### Method 1: Inline in Recipe

```yaml
_include:
  - breizorro/breizorro.yaml

cabs:
  breizorro:
    image:
      version: '0.2.1'  # Use specific version

my-recipe:
  steps:
    create_mask:
      cab: breizorro
      params:
        restored-image: test.fits
        outfile: mask.fits
```

### Method 2: Separate Override File

Create `version-override.yaml`:
```yaml
cabs:
  breizorro:
    image:
      version: updates  # Or '0.2.1', 'main', etc.
```

Run with:
```bash
stimela run -S breizorro::breizorro.yaml version-override.yaml \
  restored-image=test.fits outfile=mask.fits
```

### Method 3: Alternative Registry

Use a different container registry:
```yaml
cabs:
  breizorro:
    image:
      registry: quay.io
      name: myorg/breizorro
      version: custom
```

### Method 4: Local SIF Image

Use a locally built Singularity image:
```yaml
cabs:
  breizorro:
    image:
      path: /path/to/breizorro.sif
```

## Building Custom Images

### From Docker

1. Build Docker image:
   ```bash
   docker build -t breizorro:custom .
   ```

2. Convert to Singularity:
   ```bash
   singularity build breizorro-custom.sif docker-daemon://breizorro:custom
   ```

3. Use in recipe:
   ```yaml
   cabs:
     breizorro:
       image:
         path: /path/to/breizorro-custom.sif
   ```

### From Git Branch

Stimela can build images from git directly (if configured):
```yaml
cabs:
  breizorro:
    image:
      version: my-feature-branch
```

## Extras and Features

### Install with All Features

Default `latest` version includes `[all]` extras (catalog + GUI):
```yaml
cabs:
  breizorro:
    image:
      version: latest  # Includes [all] extras
```

### Minimal Installation

For a minimal installation (core features only):
```yaml
cabs:
  breizorro:
    image:
      version: '0.2.1'  # Specific versions are minimal by default
```

## Version Management Best Practices

### Production Workflows

Pin to specific versions for reproducibility:
```yaml
cabs:
  breizorro:
    image:
      version: '0.2.1'  # Reproducible
```

### Development/Testing

Use `latest` or development branches:
```yaml
cabs:
  breizorro:
    image:
      version: updates  # Latest features
```

### Multi-Version Testing

Test with multiple versions in CI:
```yaml
# test-v0.2.1.yaml
cabs:
  breizorro:
    image:
      version: '0.2.1'

# test-latest.yaml  
cabs:
  breizorro:
    image:
      version: latest
```

Run tests:
```bash
stimela run -S recipe.yml test-v0.2.1.yaml
stimela run -S recipe.yml test-latest.yaml
```

## Troubleshooting

### Force Rebuild

Force Stimela to rebuild an image:
```bash
stimela build -r breizorro::breizorro.yaml
```

### Use Local Development Version

For active development:
```bash
# Install locally
pip install -e .[all]

# Use native backend (no container)
stimela run -N breizorro::breizorro.yaml [params...]
```

## References

- [Stimela 2 Documentation](https://stimela.readthedocs.io/)
- [Breizorro GitHub](https://github.com/ratt-ru/breizorro)
- [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
