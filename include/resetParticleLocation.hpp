#ifndef resetParticleLocation_HPP
#define resetParticleLocation_HPP

#include "treeStructure.hpp"

void resetParticleLocation(OctreeClass *tree)
{
  ContainerClass* particles;
  OctreeClass::Iterator octreeIterator(tree);
  octreeIterator.gotoTop();
  octreeIterator.getCurrentCell()->resetToInitialState();

 	tree->forEachCell([&](CellClass* cell){
			cell->resetToInitialState();});

  do{
    octreeIterator.moveDown();
    octreeIterator.gotoLeft();
    do{
      octreeIterator.getCurrentCell()->resetToInitialState();
      if (!octreeIterator.canProgressToDown())
      {
        particles = octreeIterator.getCurrentListSrc();
        particles->resetForcesAndPotential();
      }      
    } while(octreeIterator.moveRight());
    octreeIterator.gotoLeft();
  } while(octreeIterator.canProgressToDown());
}

#endif //resetParticleLocation_HPP