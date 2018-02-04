let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/KechrisLab/MSPrep
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +20 msprep.Rmd
badd +0 ~/Projects/KechrisLab/MSPrep_development/Example.R
badd +70 R/readdata.R
badd +1 R/readdata_simple.R
badd +0 term://.//3598:R\ --no-save\ --quiet
badd +0 ~/Projects/Archive/CHC/p013_dominguez_kawasaki/build.R
badd +179 ~/Projects/Archive/CHC/p013_dominguez_kawasaki/p013_dominguez_kawasaki.Rmd
badd +93 ~/Projects/Archive/CHC/p013_dominguez_kawasaki/Rmd/rq4_coronary_artery_lesions.Rmd
badd +0 R_doc
argglobal
silent! argdel *
$argadd msprep.Rmd
set stal=2
edit R/readdata_simple.R
set splitbelow splitright
wincmd _ | wincmd |
split
1wincmd k
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe '1resize ' . ((&lines * 41 + 30) / 60)
exe 'vert 1resize ' . ((&columns * 157 + 119) / 238)
exe '2resize ' . ((&lines * 41 + 30) / 60)
exe 'vert 2resize ' . ((&columns * 80 + 119) / 238)
exe '3resize ' . ((&lines * 15 + 30) / 60)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 161 - ((31 * winheight(0) + 20) / 41)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
161
normal! 07|
wincmd w
argglobal
enew
file R_doc
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
wincmd w
argglobal
if bufexists('term://.//3598:R\ --no-save\ --quiet') | buffer term://.//3598:R\ --no-save\ --quiet | else | edit term://.//3598:R\ --no-save\ --quiet | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 10015 - ((14 * winheight(0) + 7) / 15)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
10015
normal! 0
wincmd w
exe '1resize ' . ((&lines * 41 + 30) / 60)
exe 'vert 1resize ' . ((&columns * 157 + 119) / 238)
exe '2resize ' . ((&lines * 41 + 30) / 60)
exe 'vert 2resize ' . ((&columns * 80 + 119) / 238)
exe '3resize ' . ((&lines * 15 + 30) / 60)
tabedit ~/Projects/Archive/CHC/p013_dominguez_kawasaki/build.R
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 2 - ((1 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
tabedit R/readdata.R
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 2 - ((1 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
tabedit ~/Projects/KechrisLab/MSPrep_development/Example.R
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 7 - ((6 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
7
normal! 0
tabnext 1
set stal=1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
