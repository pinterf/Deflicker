/*
    Deflicker plugin for Avisynth 2.5 -  mean intensity stabilizer
  Version 0.4 August 16,2004
  (c) 2004, A.G.Balakhnin aka Fizick
  bag@hotmail.ru

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


*/
//following  includes needed
#include "windows.h"
#include "avisynth.h"

#include "stdio.h"  // for debug
#include "math.h"  // for sqrt, fabs

#include "info.h" // for info on frame
#include "stdint.h"

class Deflicker : public GenericVideoFilter {
  // Deflicker defines the name of your filter class. 
  // This name is only used internally, and does not affect the name of your filter or similar.
  // This filter extends GenericVideoFilter, which incorporates basic functionality.
  // All functions present in the filter must also be present here.


  float percent; //  influence of previuos frame mean luma (in percent) for temporal luma smoothing
  int lag;  // max distance to base frame for temporal luma smoothing (may be negative for backward time mode)
  float noise; // noise std deviation (due to motion etc)
  int scene; // threshold for new scene (mean luma difference or std. variation doubled difference)
  int lmin;       // luma min
  int lmax;       // luma max
  int border;	// exclude border at all edges (pixels)
  int info; // show info on frame
  int debug;

  // internal parameters
  int range; // = abs(lag)
  float varnoise; // noise variation = noise*noise

  int cachecapacity; // size of cache

  int* cachelist;
  int* meancache;
  int* varcache;

  char* debugbuf;
  char* messagebuf;

  int fieldbased;
  bool isYUY2;

  int opt; // 0:C, 1:SSE

public:
  // This defines that these functions are present in your class.
  // These functions must be that same as those actually implemented.
  // Since the functions are "public" they are accessible to other classes.
  // Otherwise they can only be called from functions within the class itself.

  Deflicker(PClip _child, float _percent, int _lag, float _noise, int _scene, int _lmin, int _lmax, int _border, int _info, int _debug, IScriptEnvironment* env);
  // This is the constructor. It does not return any value, and is always used, 
  //  when an instance of the class is created.
  // Since there is no code in this, this is the definition.

  ~Deflicker();
  // The is the destructor definition. This is called when the filter is destroyed.


  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
  // This is the function that AviSynth calls to get a given frame.
  // So when this functions gets called, the filter is supposed to return frame n.
};

/***************************
 * The following is the implementation
 * of the defined functions.
 ***************************/

 //Here is the actual constructor code used
Deflicker::Deflicker(PClip _child, float _percent, int _lag, float _noise, int _scene, int _lmin, int _lmax, int _border, int _info, int _debug, IScriptEnvironment* env) :
  GenericVideoFilter(_child), percent(_percent), lag(_lag), noise(_noise), scene(_scene), lmin(_lmin), lmax(_lmax), border(_border), info(_info), debug(_debug) {
  // This is the implementation of the constructor.
  // The child clip (source clip) is inherited by the GenericVideoFilter,
  //  where the following variables gets defined:
  //   PClip child;   // Contains the source clip.
  //   VideoInfo vi;  // Contains videoinfo on the source clip.

  opt = 0; // C

#ifndef X86_64
  // no x64 asm at the moment. Todo: move to intrinsics 
  if (env->GetCPUFlags() & CPUF_SSE)
  {
    opt = 1;
  }
#endif

  if(!vi.IsYUY2() && !vi.IsY() && !vi.IsYUV() && !vi.IsYUVA())
    env->ThrowError("Deflicker: input to filter must be Y, YUV or YUY2!");
  if (vi.BitsPerComponent() > 8)
    env->ThrowError("Deflicker: only 8 bit clips supported");

  isYUY2 = vi.IsYUY2();

  range = abs(lag);

  if (percent < 0 || percent >100)
    env->ThrowError("Deflicker: percent must be from 0 to 100 !");
  //	else if (percent == 0) 	range = 0;
  //	else if (percent == 100) range = 1;
  //	else if (range <= 0) range = int( 5*log(100/(100-percent)) );

  cachecapacity = range * 2 + 2; // with some reserve

  cachelist = (int*)malloc(cachecapacity * sizeof(int));
  meancache = (int*)malloc(cachecapacity * sizeof(float));
  varcache = (int*)malloc(cachecapacity * sizeof(float));

  for (int i = 0; i < cachecapacity; i++) {
    cachelist[i] = -1; // unused, free
  }

  varnoise = noise * noise;

  debugbuf = (char*)malloc(256);
  messagebuf = (char*)malloc(64);

  // avisynth frame cache
  SetCacheHints(CACHE_WINDOW, cachecapacity);

  fieldbased = vi.IsFieldBased();  // after separatefields


}

// This is where any actual destructor code used goes
Deflicker::~Deflicker() {
  // This is where you can deallocate any memory you might have used.
  free(cachelist);
  free(meancache);
  free(varcache);
  free(debugbuf);
  free(messagebuf);

}




//****************************************************************************
//	clear un-needed cache
void clear_unnecessary_cache(int* cachelist, int cachecapacity, int ndest, int range)
{
  int i;
  for (i = 0; i < cachecapacity; i++)
  {
    if (cachelist[i] > ndest + range || cachelist[i] < ndest - range)		cachelist[i] = -1;
  }
}


//****************************************************************************
// check if  of this frame is in cache, get number of it
int get_cache_number(int* cachelist, int cachecapacity, int ndest)
{
  int i;
  int found = -1;
  for (i = 0; i < cachecapacity; i++)
  {
    if (cachelist[i] == ndest)
    {
      found = i;
    }
  }
  return found;
}

//
//****************************************************************************
//
int get_free_cache_number(int* cachelist, int cachecapacity)
{
  int i;
  int found = -1;
  for (i = 0; i < cachecapacity; i++)
  {
    if (cachelist[i] == -1)
    {
      found = i;
    }
  }
  return found;
}
//
//****************************************************************************
//

const __int64 YMask = 0x00ff00ff00ff00ff;
const __int64 UVMask = 0xff00ff00ff00ff00;
const __int64 AllMask = 0xffffffffffffffff;

//
//****************************************************************************
//
_inline void SumLine_c(const BYTE* srcp, int width, int* mean, int* varline)
{
  for (int w = 0; w < width; w += 2)
  {   // use step=2 for more fast YV12 and using the same code for YUY2
    int cur = srcp[w];
    *mean += cur;
    *varline += cur * cur; // simply quadrat, sub mean later
  }
}

//
//****************************************************************************
//
#ifndef X86_64
void SumLine_isse(const BYTE* srcp, int width, int* mean, int* varline)
{
  // no inline asm for x64 in VS2019
//	for (int w = 0; w < width; w+=2) 
//	{   // use step=2 for more fast YV12 and using the same code for YUY2
//		int cur = srcp[w];
//		*mean += cur;
//		*varline += cur*cur; // simply quadrat, sub mean later
//	}
  _asm
  {
    mov		esi, srcp
    mov		ecx, width
    mov		edi, mean
    mov		ebx, varline
    movd	mm1, [edi] // will be sum of src
    pxor	mm0, mm0 // 0
    movd	mm2, [ebx] // will be sum of src*src
    xor eax, eax
    ALIGN 16
    Loop4:
    movq	mm3, [esi + eax] // get 8 byte (4 pixels)
      pand	mm3, YMask // every second byte of src
      movq	mm4, mm3
      psadbw	mm3, mm0 // sum of 4 abs(src-0)
      paddd	mm1, mm3 // cumulative sum
      pmullw	mm4, mm4 // low word of src*src
      movq	mm5, mm4
      psllw	mm5, 8 //low bytes to high
      psadbw	mm5, mm0 // sum of abs(src*src-0) low
      psrlw	mm4, 8 // high bytes to low
      psadbw	mm4, mm0 // sum of abs(src*src-0) high
      pslld	mm4, 8 // set  sum ( high) as high
      paddd	mm2, mm4 // cumulative sum of src*src low
      paddd	mm2, mm5 // cumulative sum of src*src high
      add		eax, 8
      cmp		eax, ecx
      jl		Loop4

      movd[edi], mm1 // save sum of src
      movd[ebx], mm2 // save sum of src*src

  }
}
#endif

//
//****************************************************************************
//
_inline void CorrectYUY2_c(BYTE* dstp, const BYTE* srcp, int src_width, short mult256w, short addw, short lmin, short lmax)
{ // correct luma 
  for (int w = 0; w < src_width; w += 2) // every 2
  {
    int cur = ((srcp[w] * mult256w) >> 8) + addw;
    dstp[w] = min(max(cur, lmin), lmax); //Y
    dstp[w + 1] = srcp[w + 1]; // U,V
  }
}

//
//****************************************************************************
//
#ifndef X86_64
void CorrectYUY2_isse(BYTE* dstp, const BYTE* srcp, int src_width, short mult256w, short addw, short lmin, short lmax)
{ // correct luma 
//	for (int w = 0; w < src_width; w+=2) // every 2
//	{    
//		int cur = ((srcp[w]*mult256w)>>8) + addw;  
//		dstp[w] = min(max(cur,lmin),lmax);
//	}
  _asm
  {
    mov		esi, srcp
    mov		edi, dstp
    mov		ecx, src_width
    pinsrw	mm0, mult256w, 0
    pshufw	mm0, mm0, 0 // mult256w to 4 words
    pinsrw	mm1, addw, 0
    pshufw	mm1, mm1, 0 // addw in 4 words
    pinsrw	mm2, lmin, 0
    pshufw	mm2, mm2, 0 // lmin to 4 words (packed_low)
    pinsrw	mm3, lmax, 0
    pshufw	mm3, mm3, 0 // lmax to 4 words (pached_high)
    mov		ax, 0x8000
    pinsrw	mm6, ax, 0
    pshufw	mm6, mm6, 0 //0x8000800080008000 // packed min
    movq	mm4, mm2
    paddw	mm4, mm6 // low_us
    paddw	mm3, mm6 // high_us
    pcmpeqb	mm7, mm7 // 11111111111 (packed_usmax)
    psubusw mm7, mm3 // (packed_usmax - high_us)
    paddw	mm4, mm7 // (packed_usmax - high_us + low_us)
    pxor	mm3, mm3 // 0
    xor eax, eax
    ALIGN 16
    Loop4:
    movq	mm5, [esi + eax] // get 8 bytes (4 pixels) from source
      movq	mm3, mm5
      pand	mm5, YMask
      pand	mm3, UVMask
      psllw	mm5, 8 // src*256
      pmulhuw	mm5, mm0 //  high words  of (src*256)*mult256  
      paddsw	mm5, mm1 // + addw
      // clip
      paddw	mm5, mm6 // 0x8000800080008000
      PADDUSW	MM5, mm7 //mm5+(PACKED_USMAX - HIGH_US)
      PSUBUSW MM5, mm4 //(mm5-PACKED_USMAX - HIGH_US + LOW_US)
      PADDW   MM5, mm2 // +PACKED_LOW
      por		mm5, mm3 // add UV
      movq[edi + eax], mm5 // put 8 bytes to dest
      add		eax, 4
      cmp		eax, ecx
      jl		Loop4

  }
}
#endif

//
//****************************************************************************
//
_inline void CorrectYV12_c(BYTE* dstp, const BYTE* srcp, int src_width, short mult256w, short addw, short lmin, short lmax)
{ // correct luma 
  for (int w = 0; w < src_width; w += 1) // every 1
  {
    int cur = ((srcp[w] * mult256w) >> 8) + addw;
    dstp[w] = min(max(cur, lmin), lmax);
  }
}

//
//****************************************************************************
//
#ifndef X86_64
void CorrectYV12_isse(BYTE* dstp, const BYTE* srcp, int src_width, short mult256w, short addw, short lmin, short lmax)
{ // correct luma 
//	for (int w = 0; w < src_width; w+=1) // every 1
//	{    
//		int cur = ((srcp[w]*mult256w)>>8) + addw;  
//		dstp[w] = min(max(cur,lmin),lmax);
//	}
  _asm
  {
    mov		esi, srcp
    mov		edi, dstp
    mov		ecx, src_width
    pinsrw	mm0, mult256w, 0
    pshufw	mm0, mm0, 0 // mult256w to 4 words
    pinsrw	mm1, addw, 0
    pshufw	mm1, mm1, 0 // addw in 4 words
    pinsrw	mm2, lmin, 0
    pshufw	mm2, mm2, 0 // lmin to 4 words (packed_low)
    pinsrw	mm3, lmax, 0
    pshufw	mm3, mm3, 0 // lmax to 4 words (pached_high)
    mov		ax, 0x8000
    pinsrw	mm6, ax, 0
    pshufw	mm6, mm6, 0 //0x8000800080008000 // packed min
    movq	mm4, mm2
    paddw	mm4, mm6 // low_us
    paddw	mm3, mm6 // high_us
    pcmpeqb	mm7, mm7 // 11111111111 (packed_usmax)
    psubusw mm7, mm3 // (packed_usmax - high_us)
    paddw	mm4, mm7 // (packed_usmax - high_us + low_us)
    pxor	mm3, mm3 // 0
    xor eax, eax
    ALIGN 16
    Loop4:
    movd	mm5, [esi + eax] // get 4 bytes from source
      punpcklbw mm5, mm3		// src convert to 4 words
      psllw	mm5, 8 // src*256
      pmulhuw	mm5, mm0 //  high words  of (src*256)*mult256  
      paddsw	mm5, mm1 // + addw
      // clip
      paddw	mm5, mm6 // 0x8000800080008000
      PADDUSW	MM5, mm7 //mm5+(PACKED_USMAX - HIGH_US)
      PSUBUSW MM5, mm4 //(mm5-PACKED_USMAX - HIGH_US + LOW_US)
      PADDW   MM5, mm2 // +PACKED_LOW                         
      packuswb mm5, mm5
      movd[edi + eax], mm5 // put 4 bytes to dest
      add		eax, 4
      cmp		eax, ecx
      jl		Loop4

  }
}
#endif

//
//****************************************************************************
//

static void copy_plane(PVideoFrame& destf, PVideoFrame& currf, int plane, IScriptEnvironment* env) {
  const uint8_t* srcp = currf->GetReadPtr(plane);
  int src_pitch = currf->GetPitch(plane);
  int height = currf->GetHeight(plane);
  int row_size = currf->GetRowSize(plane);
  uint8_t* destp = destf->GetWritePtr(plane);
  int dst_pitch = destf->GetPitch(plane);
  env->BitBlt(destp, dst_pitch, srcp, src_pitch, row_size, height);
}

PVideoFrame __stdcall Deflicker::GetFrame(int ndest, IScriptEnvironment* env) {
  // This is the implementation of the GetFrame function.
  // See the header definition for further info.

  PVideoFrame src;// = child->GetFrame(ncur, env);
  const BYTE* srcp;//= src->GetReadPtr();
  int src_pitch;// = src->GetPitch();
  int src_width;// = src->GetRowSize();
  int src_height;// = src->GetHeight();
  int w, h;
  int ncur, n, nbase;
  int cur;
  int mean, var, meancur, varline;
  float meansmoothed;
  float a, b, alfa, beta, var_y;
  float mult, add;
  //	int mult256i, add256i;
  int newbase;
  int lagsign; // sign of the lag

  int lastsrc = -1;

  int borderw;
  if (isYUY2)
    borderw = border * 2; // for YUY2
  else
    borderw = border; // for YV12

  if (range == 0)
  { // nothing to do, null transform ( may be  some intra-frame auto-gain in future versions)
    src = child->GetFrame(ndest, env);
    return src;
  }

  if (lag > 0) lagsign = 1;
  else if (lag < 0) lagsign = -1;  // lag <0
  else lagsign = 0; // lag=0: but it must be null transform above and return!


  clear_unnecessary_cache(cachelist, cachecapacity, ndest, range);




  // This code deals with YV12 colourspace where the Y, U and V information are
  // stored in completely separate memory areas

  // This colourspace is the most memory efficient but usually requires 3 separate loops
  // However, it can actually be easier to deal with than YUY2 depending on your filter algorithim

  // So first of all deal with the Y Plane


  nbase = ndest - range * lagsign; // base is some prev or some next
  if (nbase < 0) nbase = 0;
  if (nbase > vi.num_frames - 1) nbase = vi.num_frames - 1;

  newbase = -1;

  for (ncur = ndest; ncur * lagsign >= nbase * lagsign; ncur -= lagsign)
  {

    // get min, max, and mean luma of current frame
    // check cache
    n = get_cache_number(cachelist, cachecapacity, ncur);
    if (n >= 0)
    { // from cache
      mean = meancache[n];
      var = varcache[n];
    }
    else
    { // will calculate now
      src = child->GetFrame(ncur, env);
      lastsrc = ncur;
      // Request a Read pointer from the source frame.
      srcp = src->GetReadPtr();
      src_pitch = src->GetPitch();
      src_width = src->GetRowSize();
      src_height = src->GetHeight();


      mean = 0; // mean luma
      var = 0; // luma variation
      srcp += border * src_pitch + borderw; // start offset (skip border lines)

      if (opt == 0) {
        for (h = border; h < src_height - border; h += 1)
        {   // Loop from top line to bottom line 
          varline = 0; // line luma variation 
          SumLine_c(srcp, src_width - borderw * 2, &mean, &varline);
          var += (varline >> 5); //   /32 approximation, to not overflow
          srcp += src_pitch;   // Add the pitch of one line (in bytes) to the source image.
        }
      }
#ifndef X86_64
      else {
        for (h = border; h < src_height - border; h += 1)
        {   // Loop from top line to bottom line 
          varline = 0; // line luma variation 
//					for (w = 0; w < src_width-borderw*2; w+=2) // chanded in v.0.4
//					{   // use step=2 for more fast YV12 and using the same code for YUY2
//						cur = srcp[w];
//						mean += cur;
//						varline += cur*cur; // simply quadrat, sub mean later
//					}
          SumLine_isse(srcp, src_width - borderw * 2, &mean, &varline);
          var += (varline >> 5); //   /32 approximation, to not overflow
          srcp += src_pitch;   // Add the pitch of one line (in bytes) to the source image.
        }
        _asm emms
      }
#endif

      mean = mean / ((src_height - 2 * border) * (src_width / 2 - borderw)); // norm mean value (*2 every 2) - chanded in v.0.4
      var = var / (((src_height - 2 * border) * (src_width / 2 - borderw)) / (32)); // norm variation value (*2 every 2) - // chanded in v.0.4
//			var = var; 
      var = var - mean * mean; // correction

      // put to cache
      n = get_free_cache_number(cachelist, cachecapacity);
      cachelist[n] = ncur;
      meancache[n] = mean; // non-smoothed measured value
      varcache[n] = var;
    }
    // check scenechange

    if (var <= varnoise)
    {			// check if noise variation correct
      if (ncur != ndest)
      { //then may be next as new scene
        newbase = ncur + lagsign;// new scene at ncur+1
      }
      else {
        newbase = ncur;
      }
      break; // bad frame, use next as new scene
    }

    if (ncur != ndest)
    { //then may be next as new scene
      n = get_cache_number(cachelist, cachecapacity, ncur + lagsign); // check  frame ncur+1 (previosly calculated)
      if (abs(meancache[n] - mean) > scene || 2 * fabs(sqrt(float(varcache[n])) - sqrt(float(var))) > scene)
      {
        // we compared mean luma and std.var with previous for new scene detect
        newbase = ncur + lagsign;// new scene at ncur+1
        break;
      }
    }
  }

  if (newbase >= 0) nbase = newbase;
  // calculate initial state (for base frame)
  n = get_cache_number(cachelist, cachecapacity, nbase);

  // some intial values for model parameters:
  meansmoothed = meancache[n];//  mean luma
  var = varcache[n]; // variation
  // initial coeff.
  mult = 1;
  add = 0;
  a = 1; // multiplicative parameter
  b = 0; // additive parameter


 // we use luma stablization method from AURORA 

 // for frames from base+1 to ndest
  if ((nbase + lagsign) * lagsign <= ndest * lagsign && var > varnoise)
  { // and if base frame is not bad 

    for (ncur = nbase + lagsign; ncur * lagsign <= ndest * lagsign; ncur += lagsign)
    {
      // I use simplified AURORA method of Intensity flicker correction

      //Restoration of Archived Film and Video
      //Van Roosmalen, Peter Michael Bruce
      //Thesis Delft University of Technology - with ref. - with Summary in Dutch
      //Printed by Universal Press, Veenendaal
      //Cover photo: Philip Broos
      //ISBN 90-901-2792-5
      //Copyright © 1999 by P.M.B. van Roosmalen
      // File "Restoration of Archived Film and Video 1999.pdf"

      // Simplification: applied to whole frame, globally, without regions.

      // degradation model of process :
      // z(n) = alfa*y(n) + beta(n) + eta(n)
      // z - measured luma
      // y - true luma
      // alfa - multiplicative flicker
      // beta - additive flicker
      // eta - noise

      // solution (estimation):
      // y(n) = a*z(n) + b(n)
      // a - multiplicative parameter
      // b - additive parameter

      //assumption:
      // Var(y(n)) variation is previuos yprev variation estimation,
      // what is really var(a*zprev+b)=a*a*var(zprev)
      var_y = a * a * var; // variation of true y from old varequ and a,b

      n = get_cache_number(cachelist, cachecapacity, ncur);
      meancur = meancache[n]; //  luma partially equlized mean value
      var = varcache[n]; // new (current) variation of z
      // var eta = varnoise
      alfa = sqrt(float(var - varnoise) / var_y); // alfa estimation

     // E(z(n)) = meancur; // mean of z(n) is simply mean of luma

     //assumption:
     // E(y(n)) expectation is  percent mean y estimation:
     // E(y(n-1)) = meansmoothed;  mean of y is percent y mean  estimation

     // estimation solution:
      beta = meancur - alfa * meansmoothed;

      a = sqrt(float(var - varnoise) / (var * alfa)); // the article has error: no sqrt!
      b = -beta / alfa + meancur * (1 / alfa - a); // there was error too !

      // that is all,
      // pure solution:
      meansmoothed = a * meancur + b;

      //  but we can improve predicted solution by corrector,
      //using calculated values for iteration
      var_y = a * a * var; // variation of true y
      alfa = sqrt(float(var - varnoise) / var_y); // alfa correction
      beta = meancur - alfa * meansmoothed;  // alfa correction
      a = sqrt(float(var - varnoise) / (var * alfa)); // the article has error, without sqrt!
      b = -beta / alfa + meancur * (1 / alfa - a); // there was error too !

      // prepare  for next frame calculation

      // pure solution:
      // meansmoothed = a*meancur + b;

      // But we use some factor for stability
      // y = k*(a*z+b) + (1-k)*z
      // k = percent/100
      // good k=0.85

      // so modified stable solution for mean:
      // meansmoothed = (percent/100)*(a*meancur + b) + (1-percent/100)*meancur;

      // this may be rewrited:
      // meansmoothed = meancur*(1 + (a-1)*percent/100) + b*percent/100;

      // But the better we will use new mult and add:

      mult = 1 + (a - 1) * percent / 100;
      add = b * percent / 100;

      // final for current frame, for next iteration:
      meansmoothed = meancur * mult + add;

    }
    // we will use final summary mult and add for every original pixel

  }


  n = get_cache_number(cachelist, cachecapacity, ndest);

  if (debug != 0)
  {
    sprintf(debugbuf, "Deflicker: dest=%d base=%d mean=%d sq.var=%5.1f mult=%5.3f add=%6.1f smoothed=%5.1f\n", ndest, nbase, meancache[n], sqrt(float(varcache[n])), mult, add, meansmoothed);
    OutputDebugString(debugbuf);
  }

  // get short integer scaled coefficients
//	mult256i = mult*256;
//	add256i = add*256;
  short mult256w = mult * 256;
  short addw = add;


  // now make correction
  if (lastsrc != ndest)
  {
    src = child->GetFrame(ndest, env); // get frame pointer if was not last processed
  }
  // Construct a frame based on the information of the current frame
  // contained in the "vi" struct.
  PVideoFrame dst = env->NewVideoFrame(vi);

  // Request a Read pointer from the source frame.
  srcp = src->GetReadPtr();
  src_pitch = src->GetPitch();
  src_width = src->GetRowSize();
  src_height = src->GetHeight();


  // Request a Write pointer from the newly created destination image.
  BYTE* dstp = dst->GetWritePtr();

  // Requests pitch (length of a line) of the destination image.
  const int dst_pitch = dst->GetPitch();

  // Requests rowsize (number of used bytes in a line.
  const int dst_width = dst->GetRowSize();

  // Requests the height of the destination image.
  const int dst_height = dst->GetHeight();

  if (!isYUY2)
  { 
    // copy chroma planes
    if (vi.NumComponents() >= 3) {
      copy_plane(dst, src, PLANAR_U, env);
      copy_plane(dst, src, PLANAR_V, env);
    }
    // copy alpha plane
    if(vi.NumComponents() == 4)
      copy_plane(dst, src, PLANAR_A, env);
  }

  // deal with the Y Plane

  if (opt == 0) { // C
    // luma correction	
    for (h = 0; h < src_height; h++)
    {       // Loop from top line to bottom line 
      if (isYUY2)
      {
        CorrectYUY2_c(dstp, srcp, src_width, mult256w, addw, lmin, lmax);
      }
      else
      { // YV12
        CorrectYV12_c(dstp, srcp, src_width, mult256w, addw, lmin, lmax);
      }
      srcp += src_pitch;            // Add the pitch  of one line (in bytes) to the source image.
      dstp += dst_pitch;            // Add the pitch of one line (in bytes) to the dest image.
    }
  }
#ifndef X86_64
  else {

    // luma correction	
    for (h = 0; h < src_height; h++)
    {       // Loop from top line to bottom line 
      if (isYUY2)
      {
        //			for (w = 0; w < src_width; w+=2) // every 2
        //			{    
        //				cur = (srcp[w]*mult256w)>>8 + addw;  // scaled by 256 integer - fast 
        //				dstp[w] = min(max(cur,lmin),lmax);
        //				dstp[w+1] = srcp[w+1];
        //			}
        CorrectYUY2_isse(dstp, srcp, src_width, mult256w, addw, lmin, lmax);
      }
      else
      { // YV12
//			for (w = 0; w < src_width; w++) // all
//			{    
//				cur = (srcp[w]*mult256w )>>8 + addw;  // scaled by 256 integer - fast 
//				dstp[w] = min(max(cur,lmin),lmax);
//			}
        CorrectYV12_isse(dstp, srcp, src_width, mult256w, addw, lmin, lmax);
      }
      srcp += src_pitch;            // Add the pitch  of one line (in bytes) to the source image.
      dstp += dst_pitch;            // Add the pitch of one line (in bytes) to the dest image.
    }
    __asm emms
  }
#endif

  // make border visible
  if (info != 0 && border > 0)
  {// show border
    dstp = dst->GetWritePtr();
    // top
    dstp += border * dst_pitch;
    if (isYUY2)
    {
      for (w = borderw; w < src_width - borderw; w += 2)
      {
        dstp[w] = lmax + lmin - dstp[w];
      }
      dstp += dst_pitch;
      // left, right
      for (h = border + 1; h < src_height - border - 1; h++)
      {
        dstp[borderw] = lmax + lmin - dstp[borderw]; // left border
        dstp[src_width - borderw - 2] = lmax + lmin - dstp[src_width - borderw - 2]; //right border
        dstp += dst_pitch;
      }
      // bottom border
      for (w = borderw; w < src_width - borderw; w += 2)
      {
        dstp[w] = lmax + lmin - dstp[w];
      }
    }
    else
    { // YV12
      for (w = border; w < src_width - border; w++)
      {
        dstp[w] = lmax + lmin - dstp[w];
      }
      dstp += dst_pitch;
      // left, right
      for (h = border + 1; h < src_height - border - 1; h++)
      {
        dstp[border] = lmax + lmin - dstp[border]; // left border
        dstp[src_width - border - 1] = lmax + lmin - dstp[src_width - border - 1]; //right border
        dstp += dst_pitch;
      }
      // bottom border
      for (w = border; w < src_width - border; w++)
      {
        dstp[w] = lmax + lmin - dstp[w];
      }
    }
  }

  int xmsg = (fieldbased != 0) ? (dst_width / 4 - 15 + (ndest % 2) * 15) : dst_width / 4 - 8; // x-position of odd fields message
  int ymsg = dst_height / 40 - 4;  // y position
  if (info != 0)
  { // show text info on frame
    // for YUY2 the text is not very good, but quite readable :-)
    if (nbase != ndest)
    {
      sprintf(messagebuf, " DeFlicker");
      DrawString(dst, vi, xmsg, ymsg, messagebuf);
    }
    else
    {
      sprintf(messagebuf, " DeFlicker:NEW"); // new scene
      DrawString(dst, vi, xmsg, ymsg, messagebuf);
    }
    sprintf(messagebuf, " frame=%7d", ndest);
    DrawString(dst, vi, xmsg, ymsg + 1, messagebuf);
    sprintf(messagebuf, " base =%7d", nbase);
    DrawString(dst, vi, xmsg, ymsg + 2, messagebuf);
    sprintf(messagebuf, " mean   =%3d", meancache[n]);
    DrawString(dst, vi, xmsg, ymsg + 3, messagebuf);
    sprintf(messagebuf, " sq.var =%5.1f", sqrt(float(varcache[n])));
    DrawString(dst, vi, xmsg, ymsg + 4, messagebuf);
    sprintf(messagebuf, " mult= %7.3f", mult);
    DrawString(dst, vi, xmsg, ymsg + 5, messagebuf);
    sprintf(messagebuf, " add = %7.1f", add);
    DrawString(dst, vi, xmsg, ymsg + 6, messagebuf);
    sprintf(messagebuf, " m.fixed=%5.1f", meansmoothed);
    DrawString(dst, vi, xmsg, ymsg + 7, messagebuf);
  }


  // end of Y plane Code

  // As we now are finished processing the image, we return the destination image.
  return dst;
}


// This is the function that created the filter, when the filter has been called.
// This can be used for simple parameter checking, so it is possible to create different filters,
// based on the arguments recieved.

AVSValue __cdecl Create_Deflicker(AVSValue args, void* user_data, IScriptEnvironment* env) {
  return new Deflicker(args[0].AsClip(), // the 0th parameter is the source clip
    args[1].AsFloat(85), // percent -previuos frame influence factor (percent) 
    args[2].AsInt(25), //lag - max distance to base frame
    args[3].AsFloat(10), // noise - noise std deviation 
    args[4].AsInt(40), //scene - threshold for new scene (luma levels)
    args[5].AsInt(0), //lmin - luma min
    args[6].AsInt(255), //lmax -  luma max
    args[7].AsInt(0), //border - exclude from estimation the border at all edges (pixels)
    args[8].AsBool(false), //info
    args[9].AsBool(false), //debug
    env);
  // Calls the constructor with the arguments provied.
}


// The following function is the function that actually registers the filter in AviSynth
// It is called automatically, when the plugin is loaded to see which functions this filter contains.

const AVS_Linkage* AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors) {
  AVS_linkage = vectors;
  env->AddFunction("Deflicker", "c[percent]f[lag]i[noise]f[scene]i[lmin]i[lmax]i[border]i[info]b[debug]b", Create_Deflicker, 0);
  return "Deflicker plugin";
}

