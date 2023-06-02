# ztfrest
ZTF REaltime Search and Triggering 

## Reference
 
[Andreoni & Coughlin et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...918...63A/abstract); The Astrophysical Journal, Volume 918, Issue 2, id.63, 16 pp.

## Project outline

In this project we query the ZTF Alert database using Kowalski and we use forced photometry to identify rapidly fading candidates in quasi real time and trigger follow-up obseervations.

![Transient searching flowchart](flowchart_paper.png)

## Software
* Python, astropy, pandas, matplotlib, Slack
* Kowalski [(Duev et al., 2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.3582D/abstract) - GitHub: https://github.com/dmitryduev/kowalski
* ForcePhotZTF [(Yao et al., 2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...886..152Y/abstract) - GitHub: https://github.com/yaoyuhan/ForcePhotZTF

