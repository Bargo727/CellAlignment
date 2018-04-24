/////////////////////////////////////////////////////////////////////////////////////////
//
// gro is protected by the UW OPEN SOURCE LICENSE, which is summarized here.
// Please see the file LICENSE.txt for the complete license.
//
// THE SOFTWARE (AS DEFINED BELOW) AND HARDWARE DESIGNS (AS DEFINED BELOW) IS PROVIDED
// UNDER THE TERMS OF THIS OPEN SOURCE LICENSE (“LICENSE”).  THE SOFTWARE IS PROTECTED
// BY COPYRIGHT AND/OR OTHER APPLICABLE LAW.  ANY USE OF THIS SOFTWARE OTHER THAN AS
// AUTHORIZED UNDER THIS LICENSE OR COPYRIGHT LAW IS PROHIBITED.
//
// BY EXERCISING ANY RIGHTS TO THE SOFTWARE AND/OR HARDWARE PROVIDED HERE, YOU ACCEPT AND
// AGREE TO BE BOUND BY THE TERMS OF THIS LICENSE.  TO THE EXTENT THIS LICENSE MAY BE
// CONSIDERED A CONTRACT, THE UNIVERSITY OF WASHINGTON (“UW”) GRANTS YOU THE RIGHTS
// CONTAINED HERE IN CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND CONDITIONS.
//
// TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION
//
//

#include "jEColi.h"
#include "Programs.h"

#define FMULT 0.125

extern double          pAlpha;
extern double          pBeta;
extern double          pDil;
extern int            g_inducerFlag;

void jEColi::compute_parameter_derivatives ( void ) {

  lambda = sqrt ( 10 * get_param ( "ecoli_growth_rate" ) * get_param ( "ecoli_growth_rate" ) \
                  / get_param ( "ecoli_division_size_variance" ) );

//jMODS:
  if(0)
//  if(1)
//  if( (1 == g_inducerFlag)  )
//  if( (1 == g_inducerFlag) && (ECOLITYPE_02 == cellType) )
          div_vol = 2.0;//2/3 ish size
      else
          div_vol = 3.14;//=4um*pi*(1/2um)^2...thus, divide at length = 4um

//  div_vol = 3.5;
//  div_vol = get_param ( "ecoli_division_size_mean" ) - get_param ( "ecoli_growth_rate" ) * 10 / lambda;
}

jEColi::jEColi (jChipmunk_Habitat *h, World * w,
              float x, float y, float a, float v ,
                jEColi_t type)
    : Cell ( w ), volume ( v )
{

//jMODS:
  parent = w;
  habitat = h;
  cellType = type;

  compute_parameter_derivatives();

  jChipmunk_EColiInit_t cell;
      cell.habitat = habitat;
      if(ECOLITYPE_STUCKCELL == cellType)
      {
          cell.mass = 1.0e6;  cell.moment = 1.0e6;
//          cell.length = 80.0;  cell.width = 50.0;
          cell.length = 80.0;  cell.width = 20.0;
      }
      else
      {
          cell.mass = 1.0;  cell.moment = 100.0;
//          cell.mass = 0.5;  cell.moment = 100.0;
//          cell.mass = 0.1;  cell.moment = 50.0;
//          cell.mass = 1.0;  cell.moment = 1.0;
          cell.length = jlengthFromVol(v);  cell.width = 10.0;
      }
      cell.x = x;  cell.y = y;  cell.angle = a;
      cell.vx = 0.0;  cell.vy = 0.0;
//30Mar.2016:  added parent pointer:
      cell.owner = this;

  cpCell = new jChipmunk_EColi(&cell);

  int i;
  for ( i=0; i<MAX_STATE_NUM; i++ ) q[i] = 0;
  for ( i=0; i<MAX_REP_NUM; i++ ) rep[i] = 0;

  div_count = 0;
  force_div = false;

  protein = 0.0;
  plasmidRate = 20;
//  plasmidRate = 10;
  plasmidCount = plasmidRate;

//  gen.seed(time(NULL));

}

jEColi::~jEColi()
{
    if ( cpCell != NULL ) {
      delete cpCell;
    }
    if ( gro_program != NULL ) {
      delete gro_program;
    }
}

Value * jEColi::eval ( Expr * e ) {

  if ( program != NULL )
    program->eval ( world, this, e );
}

#ifndef NOGUI
void jEColi::render ( Theme * theme, GroPainter * painter ) {

//jMODS:  NB:  Ecoli cell model should not have knowledge/access of its chipmunk implementation/model


//jMODS:
//    cpPolyShape * poly = (cpPolyShape *) shape;
//  int count = poly->numVerts;
  int count = cpPolyShapeGetCount(shape);

  //jMODS:  a fortiori, Ecoli cell model should not control how it is rendered (display)--it should call a cell view controller
  //  thus, this code should be a callback from the cell model.

  double
    gfp = ( rep[GFP] / volume - world->get_param ( "gfp_saturation_min" ) ) / ( world->get_param ( "gfp_saturation_max" ) - world->get_param ( "gfp_saturation_min" ) ),
    rfp = ( rep[RFP] / volume - world->get_param ( "rfp_saturation_min" ) ) / ( world->get_param ( "rfp_saturation_max" ) - world->get_param ( "rfp_saturation_min" ) ),
    yfp = ( rep[YFP] / volume - world->get_param ( "yfp_saturation_min" ) ) / ( world->get_param ( "yfp_saturation_max" ) - world->get_param ( "yfp_saturation_min" ) ),
    cfp = ( rep[CFP] / volume - world->get_param ( "cfp_saturation_min" ) ) / ( world->get_param ( "cfp_saturation_max" ) - world->get_param ( "cfp_saturation_min" ) );

  theme->apply_ecoli_edge_color ( painter, is_selected() );

  QColor col;

  col.setRgbF( qMin(1.0,rfp + yfp),
               qMin(1.0,gfp + yfp + cfp),
               qMin(1.0,cfp),
               0.75);

  painter->setBrush(col);

  QPainterPath path;


//jMODS:
#define jcomputex (body->p.x + cosa*v.x - sina*v.y)
#define jcomputey (body->p.y + sina*v.x + cosa*v.y)

  float a = body->a;
  float cosa = cos(a);
  float sina = sin(a);

//old chipmunk version:
//  cpVect v = poly->tVerts[0];
//  path.moveTo(v.x, v.y);

//  cpVect v = cpPolyShapeGetVert(shape, 0);
//  path.moveTo(jcomputex, jcomputey);
  cpVect v = cpBodyLocalToWorld(body, cpPolyShapeGetVert(shape, 0));
  path.moveTo(v.x, v.y);

  for ( int i=1; i<count; i++ ) {
//jMODS:
//old chipmunk version:
//      v = poly->tVerts[i];
//      path.lineTo(v.x, v.y);

//      v = cpPolyShapeGetVert(shape, i);
//      path.lineTo(jcomputex, jcomputey);
      v = cpBodyLocalToWorld(body, cpPolyShapeGetVert(shape, i));
      path.lineTo(v.x, v.y);

  }

  path.closeSubpath();
  painter->drawPath(path);

}
#endif

extern double g_compressionLimit;

void jEColi::jrender ( Theme * theme, GroPainter * painter ) {

  double
    gfp = ( rep[GFP] / volume - world->get_param ( "gfp_saturation_min" ) ) / ( world->get_param ( "gfp_saturation_max" ) - world->get_param ( "gfp_saturation_min" ) ),
    rfp = ( rep[RFP] / volume - world->get_param ( "rfp_saturation_min" ) ) / ( world->get_param ( "rfp_saturation_max" ) - world->get_param ( "rfp_saturation_min" ) ),
    yfp = ( rep[YFP] / volume - world->get_param ( "yfp_saturation_min" ) ) / ( world->get_param ( "yfp_saturation_max" ) - world->get_param ( "yfp_saturation_min" ) ),
    cfp = ( rep[CFP] / volume - world->get_param ( "cfp_saturation_min" ) ) / ( world->get_param ( "cfp_saturation_max" ) - world->get_param ( "cfp_saturation_min" ) );
  theme->apply_ecoli_edge_color ( painter, is_selected() );

  QColor col;

  //JW:  this is bad design to access the cpCell stuff here
  // keeping it for now...
  cpFloat compression = cpCell->compression;
  cpFloat climit = g_compressionLimit;
  float forceScale = 1.0;


  //hard-code 20 pixels = 2um is smallest length
  //maxC is thus the max. dilution #

  //  double maxC = maxP/30.0;
  //  double maxC = maxP/20.0;

  // maxP = (l/20)*a/b

    double dilution = protein/jget_length();
  //  forceScale = dilution/maxP;  //palpha is scaled by length
  //  forceScale = protein/maxP;  //palpha is scaled by length

  // = p/l *20/(a/b) = p/(a/b) * 20/l
  //  forceScale = dilution/maxC;

  double growthScale = jget_length()/20.0;//this scales alpha to growth
  double maxP = (growthScale*pAlpha)/(pDil+pBeta);//max protein at current length
//  double maxP = (growthScale*pAlpha/pBeta);//max protein at current length

  forceScale = 0.5 * protein/maxP;

//    forceScale *= 5.0;//scale dt:  FIX to sync
//  forceScale *= 20.0;//scale dt:  FIX to sync
//  forceScale *= 15.0;//scale dt:  FIX to sync

//    forceScale *= compression/climit;

  bool clipScale = false;

  if(forceScale > 1.0)
  {
      forceScale = 1.0;
      clipScale = true;
  }
  if(forceScale < 0.0)
  {
      forceScale = 0.0;
  }

   //over-ride:
//  forceScale = 1.0;

  //#define COLORSHADECELLS

#ifdef COLORSHADECELLS
{
  if (compression < climit)
  {
      forceScale = 1.0;
  }
  else
  {
      cpFloat c = (compression-climit)/climit;
      if (c > 1.0) c=1.0;
      forceScale = (1.0-c);
  }

  col.setRgbF(0.0,forceScale,0.0,0.75);
}
#endif
  //need to make a parameter from cp model:
//  float forceScale = 1.0 - cpCell->compression/(2.0*g_compressionLimit);
//  float forceScale = 1.0;//-(cpCell->springRestLengthVar)/2000;
//10 MAR. 2016:
  switch(cellType)
  {
  case bECOLI_MORAN_A:
      col.setRgbF(0.0,1.0,0.0,0.75);
      break;
  case bECOLI_MORAN_B:
      col.setRgbF(0.0,0.0,1.0,0.75);
      break;



  case ECOLITYPE_01:
//      if(1.0 == forceScale)
//        if(0)
//        if(1 == g_inducerFlag)
      if(clipScale)
          //shade less red:
//          col.setRgbF(0.75*forceScale,forceScale,0.0,0.75);
          col.setRgbF(0.75,1.0,0.0,0.75);
      else
          //YELLOW:
          col.setRgbF(forceScale,forceScale,0.0,0.75);

      //      col.setRgbF(0.0,forceScale,0.0,0.75);
//      col.setRgbF(0.0,1.0,0.0,0.75);
      //      col.setRgbF(1.0,1.0,0.0,0.75);
      //      col.setRgbF(0.75,1.0,0.0,0.75);
      break;
      case ECOLITYPE_02:
//      col.setRgbF(0.0,forceScale,forceScale,0.75);
//  CYAN:
//      col.setRgbF(0.0,1.0,1.0,0.75);
//  BLUEISH:
//      col.setRgbF(0.0,0.25,1.0,0.75);
//      if(1.0 == forceScale)
//        if(1 == g_inducerFlag)
        if(0)
          //shade cyan:
          col.setRgbF(0.0,1.0,1.0,0.75);
      else
          col.setRgbF(0.0,0.25*forceScale,forceScale,0.75);
      break;
  }

#ifdef COLORSHADECELLS
  //set cell to RED if under full shutdown
  if(0.0 == forceScale)
  {
      col.setRgbF(1.0,0.0,0.0,0.75);
  }
#endif
//  col.setRgbF( qMin(1.0,rfp + yfp),
//               qMin(1.0,gfp + yfp + cfp),
//               qMin(1.0,cfp),
//               0.75);

  QPainterPath path;

  float dx = cpCell->shapeLength/2;
  float dy = (cpCell->shapeHeight/2)*1.5;
  cpVect v;

//replaced with code with cell ends:
#define DEBUG_DRAWING false

if(DEBUG_DRAWING)
{
    QColor col2;
    col2.setNamedColor("white");
  painter->setBrush(col2);

  //upper-left (-y is up)
  v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx, -dy));
  path.moveTo(v.x, v.y);
  //lower-left
      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx, dy));
      path.lineTo(v.x, v.y);

//  /* if draw half-cells, uncomment this:
      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(dx, dy));
      path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(dx, -dy));
      path.lineTo(v.x, v.y);

      v = cpBodyLocalToWorld(cpCell->bodyB, cpv(-dx, -dy));
      path.moveTo(v.x, v.y);
          v = cpBodyLocalToWorld(cpCell->bodyB, cpv(-dx, dy));
          path.lineTo(v.x, v.y);
//  to here */
      //lower-right
          v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx, dy));
          path.lineTo(v.x, v.y);
      //upper-right
          v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx, -dy));
          path.lineTo(v.x, v.y);

  path.closeSubpath();
  painter->drawPath(path);
}//end if(0)

    painter->setBrush(col);

//21Mar.2016:  draw cell ends:
  float dx1 = dx + cpCell->shapeHeight/3.0;
  float dx2 = dx + dy;
  float dy1 = dy - cpCell->shapeHeight/6.0;
  float dy2 = cpCell->shapeHeight/6.0;

//LEFT HALF OF CELL:
  //upper-left (-y is up)
  v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx, -dy));
  path.moveTo(v.x, v.y);

      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx1, -dy1));
  path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx2, -dy2));
  path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx2, dy2));
  path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx1, dy1));
  path.lineTo(v.x, v.y);

  //lower-left
    v = cpBodyLocalToWorld(cpCell->bodyA, cpv(-dx, dy));
  path.lineTo(v.x, v.y);

//RIGHT HALF OF CELL:
  //lower-right
    v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx, dy));
  path.lineTo(v.x, v.y);

      v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx1, dy1));
  path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx2, dy2));
  path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx2, -dy2));
  path.lineTo(v.x, v.y);
      v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx1, -dy1));
  path.lineTo(v.x, v.y);

  //upper-right
      v = cpBodyLocalToWorld(cpCell->bodyB, cpv(dx, -dy));
  path.lineTo(v.x, v.y);

  path.closeSubpath();
  painter->drawPath(path);

}

 int jEColi::jupdate( void ) {

     //call the chipmunk model--pass growth rate as a parameter:
     int status;
     status = cpCell->updateModel(1.0);

//     if(1 == g_inducerFlag)
//         cpCell->applyScaledCompressionForce(0.06);


    if(-1 == status)
     {
        mark_for_death();
        return -1;
     }

     volume = jget_volume();
     //volume += get_param ( "ecoli_growth_rate" ) * rand_exponential ( 1 / world->get_sim_dt() ) * volume;

     if ( volume > div_vol )
//         if ( volume > div_vol && frand() < lambda * world->get_sim_dt() )
        div_count++;

  return 0;
}

jEColi * jEColi::jdivide ( void ) {

//jMODS:  not sure about the div count need:
//    if ( div_count >= 10 || force_div ) {
//    if (0 ) {
    if ( div_count >= 1 || force_div ) {

    int r = frand() > 0.5 ? 1 : -1;

    div_count = 0;
    force_div = false;

    //equal distribution +/- 0.05
    float frac = 0.5 + 0.1 * ( frand() - 0.5 );
    float oldvol = volume;

    double oldProtein = protein;
    protein = frac*oldProtein;

    int oldPlasmid = plasmidCount;
    int splitPlasmid = 0;

    //binomial distribution of plasmids:
    for(int i=0; i<plasmidCount; i++)
    {
        if (0.5 > frand())
            splitPlasmid++;
    }
    //update this daughter cell:
//    plasmidCount = splitPlasmid;
    //use even splitting of plasmid, with random odd-plasmid daugther
    plasmidCount = oldPlasmid/2;
    if (0.5 > frand()) plasmidCount = oldPlasmid - plasmidCount;
 //============================================================
    //this should all move to the cp model:
    // like cpModel_divideCell(); --with perhaps the rand. parameters

    jChipmunk_EColiInit_t cell;

    volume = frac * oldvol;
    float a = cpCell->angle;
//    float da = 0.5 * (frand()-0.5);
//    float da = 0.20 * (frand()-0.5);
    float da = 0.0;

    float oldsize = cpCell->length;
    cpVect oldpos = cpCell->center;

    cpVect newpos = cpCell->center + cpvmult ( cpv ( cos ( a - r*da ), sin ( a - r*da ) ),
                                               (-r)*0.5*oldsize*(1-frac) );
//    cpVect vel = cpvzero;
    cpVect vel = (1 == r) ?  cpCell->velB:cpCell->velA;

    cell.habitat = habitat;
    //eventually read these from the parent cell:
    cell.mass = cpCell->initialMass;
    cell.moment =cpCell->initialMoment;

    cell.x = newpos.x;  cell.y = newpos.y;  cell.angle = a - r*da;
    cell.length = jlengthFromVol(volume);  cell.width = 10.0;
//    cell.vx = 0.0;  cell.vy = 0.0;
    cell.vx = vel.x;  cell.vy = vel.y;
//    cell.vx = cpCell->velA.x;  cell.vy = cpCell->velA.y;
    cell.owner = this;
    jChipmunk_EColi *newCell = new jChipmunk_EColi(&cell);

    //set rest length: NEED TO MAKE SETTER WHICH CALLS cp FUNCTION!
    //need to use averaged values to avoid noise in daughter cells:

    float newRestLength = (cpCell->springRestLengthVar - cpCell->length) + cell.length;

//    float newRestLength = 1.0*(cpCell->springRestLengthVar - cpCell->separationDistance);
//    float newRestLength = (cpCell->springRestLengthVar);
//    float newRestLength = (cpCell->shapeLength + cpCell->dRL);
//      float newRestLength = (cpCell->shapeLength * 2.0f);

    newCell->springRestLengthVar = newRestLength;
    newCell->compression = newCell->springRestLengthVar - newCell->length;

    //determine whether it is the "right" or "left" cell;  r=1 is "right"=bodyB
    //need to use averaged valued!
    //need to split the left/right halves of each new cell also so inner half is like c.o.m. of old cell.
    cpVect newV1, newV2;
    if(1 == r)
    {
        newV1 = cpCell->velB; newV2 = cpCell->velA;
    }
    else
    {
        newV1 = cpCell->velA; newV2 = cpCell->velB;
    }
//    if( ( cpvlength(cpCell->velA) > 20.0)
//            || ( cpvlength(cpCell->velB) > 20.0))
//    {
//        newV1 = newV2 = cpvzero;
//    }
    newCell->setVelocity(newV1);

    //replace the old cell with the new one:
    delete cpCell;
    cpCell = newCell;

//============================================================
//  new daughter cell
//============================================================
    float dvol = (1-frac)*oldvol;

    //make cell #2 ("daughter") appear:
    jEColi *daughter = new jEColi (habitat, world,
                                   oldpos.x + r*0.5*oldsize*frac*cos ( a + r*da ),
                                   oldpos.y + r*0.5*oldsize*frac*sin ( a + r*da ),
                                   a+r*da, dvol,
                                   cellType);
    daughter->cpCell->setVelocity(newV2);
//    daughter->cpCell->setVelocity((-1 == r) ? cpCell->velB:cpCell->velA);

    //set rest length
    daughter->cpCell->springRestLengthVar = newRestLength - cell.length + jlengthFromVol(dvol);
    daughter->cpCell->compression = daughter->cpCell->springRestLengthVar - daughter->cpCell->length;

    daughter->protein = (1.0-frac)*oldProtein;
    daughter->plasmidCount = oldPlasmid - plasmidCount;

    daughter->set_param_map ( get_param_map() );
    daughter->compute_parameter_derivatives();

    if ( gro_program != NULL ) {
      daughter->set_gro_program ( split_gro_program ( gro_program, frac ) );
    }

    daughter->init ( q, rep, 1-frac );

    int i;

    for ( i=0; i<MAX_STATE_NUM; i++ ) q[i] = (int) ceil(frac*q[i]);
    for ( i=0; i<MAX_REP_NUM; i++ ) rep[i] = (int) ceil(frac*rep[i]);

    set_division_indicator(true);
    daughter->set_division_indicator(true);
    daughter->set_daughter_indicator(true);

    return daughter;

  } else return NULL;

}

void jEColi::setData(cpVect force, cpVect vel, int which)
{
    forceBuffer[which].f=force;
    forceBuffer[which].v=vel;
}
