# Netmap-evaluation
Companion repository for Netmap package

## Steps to reproduce
The project uses pixi as the main dependency manager, but one step uses R reticulate which depends on conda. Therefore you need conda installed for this step.

### Install pixi (if not already installed)
```
curl -fsSL https://pixi.sh/install.sh | bash
pixi run
```

### Sythetic benchmark
```
pixi run simulate-data
```
