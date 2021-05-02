
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

#ifndef PBRT_ACCELERATORS_OCTTREEACCEL_H
#define PBRT_ACCELERATORS_OCTTREEACCEL_H

// accelerators/octtreeaccel.h*
#include "pbrt.h"
#include "primitive.h"

namespace pbrt {

// OctTreeAccel Declarations
struct OctAccelNode;
struct BoundEdge;
class OctTreeAccel : public Aggregate {
  public:
    // OctTreeAccel Public Methods
    OctTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
                int isectCost = 80, int traversalCost = 1,
                Float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    Bounds3f WorldBound() const { return bounds; }
    ~OctTreeAccel();
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;

  private:
    // OctTreeAccel Private Methods
    void buildTree(int nodeNum, const Bounds3f &nodeBounds, int *primNums,
                   int nPrimitives, int depth);

    // OctTreeAccel Private Data
    const int isectCost, traversalCost, maxPrims;
    const Float emptyBonus;
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::vector<int> primitiveIndices;
    OctAccelNode *nodes;
    int nAllocedNodes, nextFreeNode;
    Bounds3f bounds;

    // CGRA408 code
    //---//
    std::vector<std::string> octLookup = {"000", "001", "010", "011",
                                          "100", "101", "110", "111"};
    std::vector<std::vector<int>> prims;

    std::vector<std::vector<int>> getPrimsAfterSplit(int *prims, const Bounds3f &nodeBounds);
    unsigned int getOctLookupValue(std::string key);
    void subdivideBound(const Bounds3f &bound, std::vector<Bounds3f> &bounds,
                        unsigned depth);
    //---//
};

struct OctToDo {
    const OctAccelNode *node;
    Float tMin, tMax;
};

std::shared_ptr<OctTreeAccel> CreateOctTreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_OCTTREEACCEL_H
