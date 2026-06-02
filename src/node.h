/////////////////////////////////////////////////////////////////////
// Xiaoxun Sun & Sven Koenig @ USC 2009
// All rights reserved
/////////////////////////////////////////////////////////////////////

#ifndef MAZEH
#define MAZEH

#include "include.h"

struct gnode;
typedef struct gnode gnode;

struct snode;
typedef struct snode snode;


struct gnode // stores info needed for each graph node
{
  long long int id;
  unsigned h1;
  unsigned h2;
  unsigned long long int key;
  unsigned gmin;      // min g2 expanded so far at this node (Lex1 dominance)
  unsigned g1min;     // min g1 expanded so far at this node (Lex2 dominance)
  /*
   * version: set to the global g_bc_version the last time this node was
   * touched by bc_boastar.  Allows O(1) lazy reset instead of an
   * O(num_gnodes) loop at the start of every bc_boastar() call.
   */
  unsigned version;
  unsigned long heapindex;
  snode *gopfirst;
  snode *goplast;
};


struct snode // BOA*'s search nodes
{
  int state;
  unsigned g1;
  unsigned g2;
  double key;
  unsigned long heapindex;
  snode *searchtree;
  snode *gopnext;
  snode *gopprev;
};

#endif
