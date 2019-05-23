## Deflicker ##

DeFlicker plugin can remove old film intensity flicker by temporal mean luma smoothing. 
Also it can  correct blinding of automatic gain control after flashes. 
The luma stabilizing not only improves visual impression, but can help to following noise reduction.

Versions
- 0.6 (20190523) by pinterf
  - AVX2 and SSE2 support (was: ISSE in win32, plain C in x64)
  - don't overflow over 8Mpix
  - uses all pixels for non-YUY2 colorspaces in correction calculation
    (was: only even rows were summed up in calculating internal average luma and dev)

- 0.5 (20190522) by pinterf
  - project move to Visual Studio 2019
  - Avisynth 2.6 interface, using headers from Avisynth+ project
  - Allow all 8 bit Y and YUV formats
  - x64 version (although C only)
- 0.4 (2004) by Fizick

### Links ###
https://github.com/pinterf/Deflicker
http://avisynth.org.ru/deflicker/deflicker04.zip
http://avisynth.org.ru/deflicker/deflicker.html
http://www.avisynth.nl/index.php/External_filters#Luma_Equalization