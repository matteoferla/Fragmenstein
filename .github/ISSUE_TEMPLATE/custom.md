---
name: Custom issue template
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---

When submitting a bug report, could you let me know not only what Python and Fragmenstein versions you are using, but also that of PyRosetta! Thanks:

```python
import pkg_resources, sys, platform

get_version = lambda name: pkg_resources.get_distribution(name).version

print(f'Python {sys.version}\n'+
      f'on {platform.system()} {platform.machine()}\n'+
      f'with Fragmenstein {get_version("fragmenstein")} and PyRosetta  {get_version("pyrosetta")}'
     )
```
