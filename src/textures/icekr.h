
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_TEXTURES_ICEKR_H
#define PBRT_TEXTURES_ICEKR_H

// textures/dots.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

namespace pbrt {

// IceKrTexture Declarations
template <typename T>
class IceKrTexture : public Texture<T> {
  public:
    // IceKrTexture Public Methods
    IceKrTexture(const T &value, std::unique_ptr<TextureMapping2D> mapping)
        : value(value), mapping(std::move(mapping)) {
        meansAndSdsS.push_back(Point2f(0.2, sd));
        meansAndSdsS.push_back(Point2f(0.9, sd));
        meansAndSdsS.push_back(Point2f(0.12, sd));
        meansAndSdsS.push_back(Point2f(0.01, sd));
        meansAndSdsS.push_back(Point2f(0.5, sd));

        meansAndSdsT.push_back(Point2f(0.05, sd));
        meansAndSdsT.push_back(Point2f(0.63, sd));
        meansAndSdsT.push_back(Point2f(0.46, sd));
        meansAndSdsT.push_back(Point2f(0.91, sd));
        meansAndSdsT.push_back(Point2f(0.09, sd));
    }
    T Evaluate(const SurfaceInteraction &si) const {
        // CGRA408 code
        //---//
        Vector2f dstdx, dstdy;
        Point2f st = mapping->Map(si, &dstdx, &dstdy);
        Float normal_s = MultiNormDist(st[0], meansAndSdsS);
        Float normal_t = MultiNormDist(st[1], meansAndSdsT);
        return (normal_s + normal_t) / 2.0f;
        //---//
    }
   
  private:
    T value;
    std::unique_ptr<TextureMapping2D> mapping;
    // Gaussian standard deviation
    Float sd = 0.05f;

    std::vector<Point2f> meansAndSdsS;
    std::vector<Point2f> meansAndSdsT;
};

IceKrTexture<Float> *CreateIceKrFloatTexture(const Transform &tex2world,
                                          const TextureParams &tp);
IceKrTexture<Spectrum> *CreateIceKrSpectrumTexture(
    const Transform &tex2world, const TextureParams &tp);

}  // namespace pbrt

#endif  // PBRT_TEXTURES_CONSTANT_H
