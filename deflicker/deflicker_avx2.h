#ifndef __DEFLICKER_AVX2_H__
#define __DEFLICKER_AVX2_H__

#include "avisynth.h"
#include "stdint.h"

void SumFrame_YUY2_avx2(const BYTE* srcp, int src_pitch, int width, int height, int borderw, int border, int64_t* mean64, int64_t* var64);
void SumFrame_avx2(const BYTE* srcp, int src_pitch, int width, int height, int borderw, int border, int64_t* mean64, int64_t* var64);
void CorrectFrame_YUY2_avx2(BYTE* dstp, int dst_pitch, const BYTE* srcp, int src_pitch, int src_height, int src_width, short mult256w, short addw, short lmin, short lmax);
void CorrectFrame_avx2(BYTE* dstp, int dst_pitch, const BYTE* srcp, int src_pitch, int src_height, int src_width, short mult256w, short addw, short lmin, short lmax);

#endif	// __DEFLICKER_AVX2_H__
