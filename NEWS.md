## Version 3.3.3

### Bug Fixes
Fixed a WARNING appearing on clang version 15.0.6, Fedora Linux 36.

---

## Version 3.3.2

### Bug Fixes 
Fixed a major bug that severely affect the computational complexity of the procedure. Should you have an older version installed (lower than 3.3.2) please make sure you update your DeCAFS package either through CRAN or GitHub. You can check your version number at the bottom of the documentation page of DeCAFS, via `help("DeCAFS")`. 

---

## Version 3.3.1

### What's new
In addition to the automatic model selection, we introduced a graphical iterative model selection procedure that aids the user in selecting an appropriate model for a given sequence of observations. This tuning procedure can seriously improve performances under more challenging scenarios. More details can be found by checking the documentation: `help("guidedModelSelection")`.