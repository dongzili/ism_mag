#git clone url
#git add code/*.py
#git add code/figure/phase_LLminusRR_mask6.5_vrange0.001_3.pdf
#git add /figure_back/demo
git add *.sh
git add code/*.py
git rm code/._*.pdf
#git rm --cached file1.txt
git commit -m "12.4 start toeplitz"
git remote add 12.4 https://github.com/fleaf5/ism_mag.git
git push -u 12.4 master
