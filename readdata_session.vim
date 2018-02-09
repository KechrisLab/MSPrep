let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/KechrisLab/MSPrep
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +0 R/readdata_simple.R
badd +265 R/readdata.R
badd +0 R/prepare.R
badd +0 term://.//88001:R\ --no-save\ --quiet
badd +0 vignette/using_MSPrep.Rmd
badd +0 R/original_read_sj.R
badd +0 R/original_read.R
argglobal
silent! argdel *
$argadd R/readdata_simple.R
set stal=2
edit R/readdata_simple.R
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 146 + 214) / 428)
exe '2resize ' . ((&lines * 35 + 49) / 98)
exe 'vert 2resize ' . ((&columns * 281 + 214) / 428)
exe '3resize ' . ((&lines * 59 + 49) / 98)
exe 'vert 3resize ' . ((&columns * 281 + 214) / 428)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 134 - ((64 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
134
normal! 018|
wincmd w
argglobal
if bufexists('R/prepare.R') | buffer R/prepare.R | else | edit R/prepare.R | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 164 - ((28 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
164
normal! 022|
wincmd w
argglobal
if bufexists('term://.//88001:R\ --no-save\ --quiet') | buffer term://.//88001:R\ --no-save\ --quiet | else | edit term://.//88001:R\ --no-save\ --quiet | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 10059 - ((58 * winheight(0) + 29) / 59)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
10059
normal! 0
wincmd w
3wincmd w
exe 'vert 1resize ' . ((&columns * 146 + 214) / 428)
exe '2resize ' . ((&lines * 35 + 49) / 98)
exe 'vert 2resize ' . ((&columns * 281 + 214) / 428)
exe '3resize ' . ((&lines * 59 + 49) / 98)
exe 'vert 3resize ' . ((&columns * 281 + 214) / 428)
tabedit R/original_read_sj.R
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe '1resize ' . ((&lines * 93 + 49) / 98)
exe 'vert 1resize ' . ((&columns * 213 + 214) / 428)
exe '2resize ' . ((&lines * 93 + 49) / 98)
exe 'vert 2resize ' . ((&columns * 214 + 214) / 428)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 238 - ((92 * winheight(0) + 46) / 93)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
238
normal! 0
wincmd w
argglobal
if bufexists('R/original_read.R') | buffer R/original_read.R | else | edit R/original_read.R | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 42 - ((41 * winheight(0) + 46) / 93)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
42
normal! 031|
wincmd w
exe '1resize ' . ((&lines * 93 + 49) / 98)
exe 'vert 1resize ' . ((&columns * 213 + 214) / 428)
exe '2resize ' . ((&lines * 93 + 49) / 98)
exe 'vert 2resize ' . ((&columns * 214 + 214) / 428)
tabedit vignette/using_MSPrep.Rmd
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
let s:l = 2 - ((1 * winheight(0) + 46) / 93)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
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
