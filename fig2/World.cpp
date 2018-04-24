/////////////////////////////////////////////////////////////////////////////////////////
//
// gro is protected by the UW OPEN SOURCE LICENSE, which is summaraized here.
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
//jjw:  comipler flag for mods:  05Jan16
#define jMODS
#ifndef NOGUI
#include <GroThread.h>
#endif
#include "Micro.h"

//jjw:  add threads code:  02Feb16
#include "jthreads_cells.h"
//add code to factor chipmunk implementation:  14Feb16
#include "jchipmunk_habitat.h"
//use some globals for now, until code factored:
static jChipmunk_Habitat *habitat;
#include "jchipmunk_trap.h"
static jChipmunk_Trap *trap;

#include "jEColi.h"
static std::list<jEColi *>  jcells;

//27Mar.2016:  added file output for cell ratios:
//#define RECORD_DATA
#ifdef RECORD_DATA
    #define     THISPARAMSET "run: cg=2.0, 600x300\n"
    #define     THISFILELOG "../tifs_03_29_16_f05.txt"
    std::ofstream myfile;
#endif
unsigned int nextTimeRecord = 0;
//#define SPEED_SIZE  1024
//#define SPEED_SIZE  8
#define SPEED_SIZE  16
//#define SPEED_SIZE  32
bool notDoneSpeedRecord = true;
//for cell speed data:
std::ofstream speedFiles[SPEED_SIZE];
int dividedCell[SPEED_SIZE];
double  cellEnergy, cEmax, cEmin;

double g_compressionLimit;
double          pAlpha;
double          pBeta;
double          pDil;
double          xAlpha;

int cycle;
int triggerPeriod = 100;
int initialTriggerTime = 150;
//int initialTriggerTime = 100;
int nextTriggerTime = initialTriggerTime;

int g_inducerFlag = 0;


//SPLICED CODE FOR MORAN MODEL:
#define bM          10
#define bN          20
#define MORAN_SPACING_MICRONS           4.0f
#define MORAN_CELL_LENGTH_MICRONS       3.5f

//RATES:

double b_h = 5.0;
double b_a = 0.0;
double b_b = 0.0;
double b_kappa = 0.0;
int errorCount = 0.0;


//NOTE:  INSURE TRAP SIZE IS BIGGER THAN (M x N) x SPACING in .gro file

jEColi *cellMatrix[bM][bN];



#ifdef NOGUI
World::World ( void ) {
#else
World::World ( GroThread *ct ) : calling_thread ( ct ) {
#endif


    prog = NULL;
    chemostat_mode = false;
    next_id = 0;
    step = 0;
    print_rate = -1;
    movie_rate = -1;
    program_initialized = false;
    gro_message = "";
    stop_flag = false;
    zoom = 1.0;
    barriers = new std::list<Barrier>;

    set_sim_dt ( DEFAULT_SIM_DT );
    set_chip_dt ( DEFAULT_CHIP_DT );

    set_chip_dt(0.001);
//            set_chip_dt(0.002);
//        set_chip_dt(0.002);
//    set_chip_dt(0.005);

//    set_chip_dt(0.01);

//    set_chip_dt(0.02);
//    set_chip_dt(0.05);

//    set_chip_dt(0.1);
//    set_chip_dt(0.2);
//    set_chip_dt(0.5);

//              pAlpha = 0.1;
//              pAlpha = 1.0;
//              pAlpha = 4.0;

    //NB:  dilution rate is approx. rate = ln2/24mins = 0.03 min^-1
    pBeta = 0.001;
//    pBeta = 0.05;

    pDil = 0.03;
    pAlpha = 1e5*pDil;

//        g_compressionLimit = 1.0;
//        g_compressionLimit = 2.0;
//    g_compressionLimit = 20.0;
//    g_compressionLimit = 2.75;//65x20um
//    g_compressionLimit = 3.0;//not too bright
//    g_compressionLimit = 4.0;//setting for 85x25um trap
    g_compressionLimit = 7.5;//setting for 100x25um trap
//    g_compressionLimit = 8.0;
//    g_compressionLimit = 6.0;
//    g_compressionLimit = 16.0;
//    g_compressionLimit = 32.0;
//    g_compressionLimit = 100.0;

    // Cells
    population = new std::list<Cell *>; // deleted in ~World

    jcells.clear();

    //jMODS:
    //jcells = new std::list<Cell *>; // deleted in ~World
    space = cpSpaceNew(); // freed in ~World


    //jMODS:
    jnext_id = 0;
    chemostat_initialized = false;
    //jMODS:
    habitat = new jChipmunk_Habitat(this);
    //jMODS:
    for (int i=0;i<NUM_JTHREADS;i++)
    {
        cellThreads[i] = new jThread(this, &cellMutex, &cellCondition, 1<<i);
    }
    cellReady = 0;
}//World()
World::~World ( void ) {

    std::list<Cell *>::iterator j;

    for ( j=population->begin(); j!=population->end(); j++ ) {
        delete (*j);
    }

    unsigned int k;

    for ( k=0; k<signal_list.size(); k++ )
        delete signal_list[k];

//jMODS:
    std::list<jEColi *>::iterator jj;
    for ( jj=jcells.begin(); jj!=jcells.end(); jj++ ) {
        delete (*jj);
    }
    if(chemostat_initialized)
    {
        delete trap;
    }
    delete habitat;
    delete population;

//jMODS:
    for (int i=0;i<NUM_JTHREADS;i++)
    {
        delete cellThreads[i];
    }
#ifdef RECORD_DATA
    if (myfile.is_open())
    {
      myfile.close();
    }
#endif
    prog->destroy(this);

}//~World()
//jMODS:  moved to Cell.cpp
//void Cell::init ( const int * q0, const int * rep0, float frac ) {
void World::init ()
{
    // Time
    t = 0.0f;
    max_val = 0.0f;
// JJW: snip: moved to habitat in constructor:

    // Default parameters. These will be over-written when/if the program
    // defines them via "set". But just in case the user does not do this,
    // they are defined here.
    set_param ( "chemostat_width", 200);
    set_param ( "chemostat_height", 200);
    set_param ( "signal_area_width", 800);
    set_param ( "signal_num_divisions", 160);
    set_param ( "population_max", 1000 );
    set_param ( "signal_grid_width", 800 );
    set_param ( "signal_grid_height", 800 );
    set_param ( "signal_element_size", 5 );
    // Program
    ASSERT ( prog != NULL );
    if ( !program_initialized ) {

        prog->init(this); // if this throws an exception
        // the program will not be initialized
        // in the next line

        program_initialized = true;

    }
#ifdef RECORD_DATA
    myfile.open(THISFILELOG, std::ios::app);
    myfile << THISPARAMSET;
    myfile.close();
#endif

    //init the cell divided flag for speed data:
    for(int i=0;i<SPEED_SIZE;i++)
    {
        dividedCell[i]=0;
    }
    //jMODS:
    int rx, ry, rw, rh;
//    int rw = get_param("chemostat_width") - 10;
//    int rh = get_param("chemostat_height") - 10;
    int trapWidth = get_param("chemostat_width");
    int trapHeight = get_param("chemostat_height");

//#define MMCELLS  3
#define MMCELLS  4

//    int cellCount = 128;

    //    int cellCount = 64;
      int cellCount = 32;
//    int cellCount = 16;

    //    int cellCount = MMCELLS;
//    int cellCount = 8;
//    int cellCount = 4;
//    int cellCount = 1;
//    int cellCount = 2;


       nextTriggerTime = initialTriggerTime;

       g_inducerFlag = 0;
//       cycle = 0;
       cycle = 1;

//01SEP2016:  new overlap scanning:

     int cellArray[cellCount];

     jEColi *firstCell;
     srand (time(NULL));

//MM layout:
     int cellys[MMCELLS] = {-100,-80,-60,-40};
//     int cellys[MMCELLS] = {-40, -20, 0, 40};
//     int cellys[MMCELLS] = {-40, -20, 0};
///*
#if 0
     for(int icell = 0;icell<cellCount;icell++)
     {
         firstCell = new jEColi(habitat, this,
//                                rx, ry, (rand()%31415)/10000.0,
                                                                0, cellys[icell], 3.1415/2.0,
//                                                                0, cellys[icell],  (rand()%31415)/10000.0,
                                        DEFAULT_ECOLI_INIT_SIZE,
//                                    (icell%2) == 0 ? ECOLITYPE_01: ECOLITYPE_02);
                                     ECOLITYPE_01);
         jadd_cell(firstCell);
     }
//*/
#else
///*
///
///

     //SPLICED CODE FOR MORAN MODEL:
     int halfM=bM/2;
     int halfN=bN/2;


     double pixelSpacing = 10.0*MORAN_SPACING_MICRONS;//in pixels= 10 x micron

     for(int i = 0;i<bM;i++)
     {
         for(int j = 0;j<bN;j++)
        {

    //         if(trapWidth > 28)     rw = trapWidth - 28; else rw = 0;
    //         if(trapHeight > 10)    rh = trapHeight - 10; else rh = 20;
    //         rx = (rand() % rw) - rw/2;
    //         ry = (rand() % rh)- rh/2;

    //         cellArray[icell] = ry;

             ry = double(pixelSpacing*(i - halfM));
             rx = double(pixelSpacing*(j - halfN));



             jEColi_t cellType = ((rand()%2) == 0) ? bECOLI_MORAN_A: bECOLI_MORAN_B;
             double initAngle = ((rand()%2) == 0) ? 3.1415/2.0 : 0.0;

             firstCell = new jEColi(habitat, this,
    //                                rx, ry, (rand()%31415)/10000.0,
                                    rx, ry, initAngle,
                                            (MORAN_CELL_LENGTH_MICRONS/2.0)*DEFAULT_ECOLI_INIT_SIZE,
                                            cellType);
             jadd_cell(firstCell);

             cellMatrix[i][j] = firstCell;
         }
     }
//*/

#endif
//OLD VERSION:
//     for(int icell = 0; icell<cellCount;icell++)
//     {
//         rx = (rand() % rw) - rw/2;
//         ry = (rand() % rh)- rh/2;
//         //seed a first cell for the testing:
//         firstCell = new jEColi(habitat, this,
// //                               0, 0, 1.0*(3.1415/2.0),
// //                               rx, ry, icell*(3.1415/2.0),
//                                rx, ry, (rand()%31415)/10000.0,
//                                        DEFAULT_ECOLI_INIT_SIZE,
// //                                       (icell%2) == 0 ? ECOLITYPE_01: ECOLITYPE_02);
//                                        (icell%2) == 0 ? ECOLITYPE_01: ECOLITYPE_02);
// //                                    ECOLITYPE_SPEEDTEST);
//         jadd_cell(firstCell);
// //        firstCell->cpCell->setVelocity(cpv(rx,ry));
//     }

    //ECOLITYPE_SPEEDTEST CELL(s):
    if(0)
    {
        firstCell = new jEColi(habitat, this,
                                       0, 0, 0.0,
                                       DEFAULT_ECOLI_INIT_SIZE,
                                       ECOLITYPE_SPEEDTEST);
        jadd_cell(firstCell);
    }
    //12Jun.2016: MM test setup:
//    cellCount=0;
//    for(int icell = 0; icell<3;icell++)
//    {
//        firstCell = new jEColi(habitat, this,
//                                   0, -rh/2+8+icell*20, (3.1415/2.0),
//                                       DEFAULT_ECOLI_INIT_SIZE,
//                                       ECOLITYPE_01);
//        jadd_cell(firstCell);
//    }


//    for(int icell = 0; icell<4;icell++)
//    {
//        firstCell = new jEColi(habitat, this,
//                                       -250 + 10*icell, icell*15, 0.0,
//                                       DEFAULT_ECOLI_INIT_SIZE,
//                                       (jEColi_t)(ECOLITYPE_GAP0 + icell));
//        jadd_cell(firstCell);
//        firstCell->cpCell->setVelocity(cpv(1500.0,0.0));
//    }

//    firstCell = new jEColi(habitat, this,
//                                   0, -140, 3.14/2.0,
//                                   DEFAULT_ECOLI_INIT_SIZE,
//                                   ECOLITYPE_01);
//    jadd_cell(firstCell);

//    //17Mar.2016:  add a "stuck cell"
//    firstCell = new jEColi(habitat, this,
//                                   0, -200, CP_PI,
//                                   DEFAULT_ECOLI_INIT_SIZE,
//                                   ECOLITYPE_STUCKCELL);
//    jadd_cell(firstCell);


    //start the worker thread(s):
    for (int i=0;i<NUM_JTHREADS;i++)
    {
        cellThreads[i]->start();
    }
    //temp vars: to be dynamic, set by "program" to set trap geometry, forces
    float *corners;
    unsigned int count;
//#define TRAP_FORCE      1.0e15
#define TRAP_FORCE      1.0e13
//#define TRAP_FORCE      1000.0f
//    if ( chemostat_mode )
        if ( 1 )
    {
        int w = get_param("chemostat_width")/2,
        h = get_param("chemostat_height")/2;
        chemostat_initialized = true;
        trap = new jChipmunk_Trap(habitat,
                                  float(w), float(h),
                                  corners, count,
                                  TRAP_FORCE, TRAP_FORCE, false);
    }

    // Chemostat
    //snip
}
void World::emit_message ( std::string str, bool clear ) {

#ifndef NOGUI
    calling_thread->emit_message ( str, clear );
#else
    std::cerr << str << "\n";
#endif

}
void World::restart ( void ) {

    std::list<Cell *>::iterator j;

    for ( j=population->begin(); j!=population->end(); j++ ) {
        delete (*j);
    }

    //jMODS:
    std::list<jEColi *>::iterator jj;
    for ( jj=jcells.begin(); jj!=jcells.end(); jj++ ) {
        delete (*jj);
    }
    jnext_id = 0;

//jMODS: deprecated
//    cpSpaceFreeChildren(space);
    //this does not free all the children:  need to modify:
//    cpSpaceFree(space);
//    delete population;

    unsigned int k;
    for ( k=0; k<signal_list.size(); k++ ) {
        signal_list[k]->zero();
    }

    init();

    gro_message = "";

}
static char buf[1024];
#ifndef NOGUI
static void drawString ( GroPainter * painter, int x, int y, const char *str) {

    painter->drawText ( x, y, str );

}
void World::render ( GroPainter * painter ) {

    std::list<Cell *>::iterator j;
    int dec = 0;
    theme.apply_background ( painter );
    painter->scale(zoom,zoom);

    if ( signal_list.size() > 0 ) {
        int i,j;
        unsigned int k;
        float r, g, b;
        QColor col;
        float sum;

        for ( i=0; i<signal_list.front()->get_numx(); i++ )
        {
            for ( j=0; j<signal_list.front()->get_numy(); j++ )
            {
                r = 0; g = 0; b = 0;
                // determine overall signal level
                sum = 0.0;
                for ( k=0; k<signal_list.size(); k++ )
                {
                    sum += signal_list[k]->get_val(i,j);
                }
                if ( sum > 0.01 )
                {
                    for ( k=0; k<signal_list.size(); k++ )
                    {
                        // accumulate color
                        theme.accumulate_color ( k, signal_list[k]->get_val(i,j), &r, &g, &b );
                    }
                }
                if ( sum > 0.01  )
                {
                    col.setRgbF(r/k,g/k,b/k);
                    painter->fillRect(
                                signal_list[0]->get_minp().x+signal_list[0]->get_dw()*i,
                                signal_list[0]->get_minp().y+signal_list[0]->get_dh()*j,
                                (int) ( signal_list[0]->get_dw() ),
                                (int) ( signal_list[0]->get_dh() ),
                                col);
                }
            }
        }
        theme.apply_message_color(painter);

        // draw cross hairs for corners of signal area
        int x = get_param ( "signal_grid_width" ) / 2,
                y = get_param ( "signal_grid_height" ) / 2;

        painter->drawLine ( x-20, y,  x+20, y );
        painter->drawLine ( x, y-20,  x, y+20 );

        painter->drawLine ( -x-20, y,  -x+20, y );
        painter->drawLine ( -x, y-20,  -x, y+20 );

        painter->drawLine ( x-20, -y,  x+20, -y );
        painter->drawLine ( x, -y-20,  x, -y+20 );

        painter->drawLine ( -x-20, -y,  -x+20, -y );
        painter->drawLine ( -x, -y-20,  -x, -y+20 );

    }

    //jMODS:
    std::list<jEColi *>::iterator jj;
    for ( jj=jcells.begin(); jj!=jcells.end(); jj++ )
    {
        (*jj)->jrender ( &theme, painter );
    }

    theme.apply_chemostat_edge_color(painter);
    painter->setBrush( QBrush() );

    QPainterPath path;

    // Chemostat
//    if ( chemostat_mode )
        if ( true )
    {
        int w = get_param("chemostat_width")/2,
                h = get_param("chemostat_height")/2;

#if 0
        path.moveTo(-w,h);
        path.lineTo(-w,-h);
//open top:
        path.moveTo(w,-h);
        path.lineTo(w,h);
//closed top:
//        path.lineTo(w,-h);
//        path.lineTo(w,h);

        painter->drawPath(path);
#endif
    }
    std::list<Barrier>::iterator b;
    for ( b=barriers->begin(); b != barriers->end(); b++ )
    {
      path.moveTo((*b).x1,(*b).y1);
      path.lineTo((*b).x2,(*b).y2);
      painter->drawPath(path);
    }

    painter->reset();
    theme.apply_message_color(painter);
//jMODS:
    sprintf ( buf, "Cells: %d, Max: %d, t = %.2f min", (int) jcells.size(), (int) get_param ( "population_max" ), t );
//    sprintf ( buf, "Cells: %d, Max: %d, t = %.2f min", (int) population->size(), (int) get_param ( "population_max" ), t );
    drawString ( painter, -painter->get_size().width()/2+10, -painter->get_size().height()/2+20, buf );
    dec = 32;

    message_handler.render(painter);

}
#endif
void World::add_cell ( Cell * c ) {

    population->push_back ( c );
    c->set_id ( next_id++ );
    //cpSpaceStep(space, get_chip_dt());

}
//jMODS:
void World::jadd_cell ( jEColi * c )
{
    jcells.push_back ( c );
    c->set_id ( jnext_id++ );
}
static bool out_of_bounds ( float x, float y ) {

//    return x < -400 || x > 400 || y < -400 || y > 400;
    return x < -jTRAP_SIZE || x > jTRAP_SIZE || y < -jTRAP_SIZE || y > jTRAP_SIZE;

}
cpVect World::chemostat_flow ( float, float y, float mag ) {

//jjw: check if above or below the trap:
#ifdef jMODS
    int h = 100;
//    int h = get_param("chemostat_height")/2;
    if ( (y>h) || (y<-h))
    {
//        return cpv ( mag, 0 );
        return cpv (400*mag, 0 );
    }
    else return cpv ( 0, 0 );
#else
    if ( y > get_param("chemostat_height")/2 ) return cpv ( mag, 0 );
    else return cpv ( 0, 0 );
#endif

}
//jMODS:
void World::readyCells(unsigned int id, unsigned int onOff)
{//we are mutex locked here on entry:
    if(0 == onOff)
    {//busy:
        cellReady &= ~id;
    }
    else
    {//done:
        cellReady |= id;
        //if all threads ready, wakeup the WorldThread
        if(THREADS_READY == cellReady)
        {
            computeDoneMutex.lock();
                computeDoneCondition.wakeAll();
            computeDoneMutex.unlock();

        }
    }
}
void World::computeCells(unsigned int begin, unsigned int end)
{
    std::list<Cell *>::iterator j, pb, pe;
    pb = population->begin();
    pe = population->begin();
    std::advance(pb, begin);
    std::advance(pe, end);


    // update each cell
    for ( j=pb; j!=pe; j++ )
    {
            (*j)->update();

        // check for divisions
        Cell * d = (*j)->divide();
        if ( d != NULL )
        {
            add_cell ( d );
        }
    }



}
static cpFloat maxCellExpSpeed = 0.0;


int dummyCount;

void World::update ( void ) {

    if ( population->size() < get_param ( "population_max" ) ) {

        prog->world_update ( this );
        std::list<Cell *>::iterator j;
//jMODS:
#ifdef _OPENMP
        program something with open mp!
#endif
/*
        while(THREADS_READY != cellReady) {;}
        cellMutex.lock();
            cellReady = 0;
            cellsSize = population->size();
            cellCondition.wakeAll();
            //lock the done mutex here:
            computeDoneMutex.lock();
        cellMutex.unlock();

        //now sleep on the done condition:
            computeDoneCondition.wait(&computeDoneMutex);
        computeDoneMutex.unlock();
*/


        double bdt = 0.01;

        int im = trunc(frand() * bM);
        int in = trunc(frand() * bN);

        jEColi *pCell = cellMatrix[im][in];
        //logic

        int tmB = (im==0) ? 1:0;
        int bmB = (im==(bM-1)) ? 1:0;
        int lnB = (in==0) ? 1:0;
        int rnB = (in==(bN-1)) ? 1:0;

        int interior = ( (im != 0) &&  (im != bM-1) && (in != 0) && (in != bN-1) ) ? 1:0;

        int bAngle = (pCell->cpCell->angle == 0.0) ? 0 : 1;
        int bType = (pCell->cellType == bECOLI_MORAN_A) ? 0 : 1;

        double aa[4];
        int iFromB = bM - im;
        int iFromR = bN - in;

        aa[0] = bdt * (bAngle*b_h*(1-bmB) * exp(-b_kappa*iFromB)  + (1-bAngle)*b_h*(1-rnB)*exp(-b_kappa*iFromR) );
        aa[1] = bdt * (bAngle*b_h*(1-tmB) * exp(-b_kappa*im)  + (1-bAngle)*b_h*(1-lnB)*exp(-b_kappa*in) );

        double A[4];
        if(1 == interior)
        {
            A[0] = ( (cellMatrix[im-1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im+1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[2] = ( (cellMatrix[im][in-1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[3] = ( (cellMatrix[im][in+1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (4.0 - (A[0] + A[1] + A[2] + A[3]) );
        }
        else if( tmB && (1-lnB) && (1-rnB) )
        {
            A[0] = ( (cellMatrix[im+1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im][in-1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[2] = ( (cellMatrix[im][in+1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (3.0 - (A[0] + A[1] + A[2]));

        }
        else if( bmB && (1-lnB) && (1-rnB))
        {
            A[0] = ( (cellMatrix[im-1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im][in-1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[2] = ( (cellMatrix[im][in+1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (3.0 - (A[0] + A[1] + A[2]));

        }
        else if( lnB && (1-tmB) && (1-bmB))
        {
            A[0] = ( (cellMatrix[im+1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im-1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[2] = ( (cellMatrix[im][in+1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (3.0 - (A[0] + A[1] + A[2]));

        }
        else if( rnB && (1-tmB) && (1-bmB))
        {
            A[0] = ( (cellMatrix[im+1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im-1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[2] = ( (cellMatrix[im][in-1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (3.0 - (A[0] + A[1] + A[2]));
                                              
        }
        else if( lnB && tmB)
        {
            A[0] = ( (cellMatrix[im+1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im][in+1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (2.0 - (A[0] + A[1]));
        }
        else if( lnB && bmB)
        {
            A[0] = ( (cellMatrix[im-1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im][in+1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (2.0 - (A[0] + A[1]));
        }
        else if( rnB && tmB)
        {
            A[0] = ( (cellMatrix[im+1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im][in-1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (2.0 - (A[0] + A[1]));
        }
        else if( rnB && bmB)
        {
            A[0] = ( (cellMatrix[im-1][in])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            A[1] = ( (cellMatrix[im][in-1])->cpCell->angle == pCell->cpCell->angle) ? 1:0;
            aa[2] = b_a*bdt * (2.0 - (A[0] + A[1]));
        }

        aa[3] = 1.0 - (aa[0] + aa[1] + aa[2]);


        double rn = frand();
        int which;

        double cumsum[4];
        for(int i=0; i<4; i++)
        {
            if(i==0)
            {
                cumsum[i] = aa[i];
            }
            else
            {
                cumsum[i] = aa[i] + cumsum[i-1];
            }
        }
        for(int i=0; i<4; i++)
        {
            cumsum[i] /= cumsum[3];
            if(cumsum[i] >= rn)
                {
                    which = i;
                    break;
                }
        }

        double newAngle;
        jEColi_t newType;

        if(2 == which)
        {
            newAngle = (bAngle == 0) ?  (3.1415/2.0) : 0.0;
            newType = (bType==0) ? bECOLI_MORAN_A:bECOLI_MORAN_B;
        }
        else if(0 == which)
        {
            if(0 == bAngle)
            {
                //for (int i = bN-1; i>in; i--)
                //{
                  //  double oldAngle = (cellMatrix[im][i-1])->cpCell->angle;
                   // cellMatrix[im][i]->cpCell->setAngle(oldAngle);
                    //cellMatrix[im][i]->cellType = cellMatrix[im][i-1]->cellType;

                    //double oldAngle = (cellMatrix[im][i-1])->cpCell->angle;
                    //cellMatrix[im][i]->cpCell->angle = cellMatrix[im][i-1]->cpCell->angle;
                    //cellMatrix[im][i]->cellType = cellMatrix[im][i-1]->cellType;
                //}
                cellMatrix[im][in+1]->cpCell->angle = cellMatrix[im][in]->cpCell->angle;
                cellMatrix[im][in+1]->cellType = cellMatrix[im][in]->cellType;
    
            }
            else if(1 == bAngle)
            {
                //for (int i = bM-1; i>im; i--)
                //{
                  //  double oldAngle = (cellMatrix[i-1][in])->cpCell->angle;
                   // cellMatrix[i][in]->cpCell->setAngle(oldAngle);
                    //cellMatrix[i][in]->cellType = cellMatrix[i-1][in]->cellType;

                    //double oldAngle = (cellMatrix[i-1][in])->cpCell->angle;
                    //cellMatrix[i][in]->cpCell->angle = cellMatrix[i-1][in]->cpCell->angle;
                    //cellMatrix[i][in]->cellType = cellMatrix[i-1][in]->cellType;
                //}
                cellMatrix[im+1][in]->cpCell->angle = cellMatrix[im][in]->cpCell->angle;
                cellMatrix[im+1][in]->cellType = cellMatrix[im][in]->cellType;
 
            }
            else
            {
                errorCount++;
            }
            newAngle = (bAngle == 0) ?  0.0 : (3.1415/2.0);
            newType = (bType == 0) ? bECOLI_MORAN_A:bECOLI_MORAN_B;
        }
        else if(1 == which)
        {
            if(bAngle == 0)
            {
               // for (int i = 0; i<in; i++)
                //{
                   // double oldAngle = (cellMatrix[im][i+1])->cpCell->angle;
                   // cellMatrix[im][i]->cpCell->setAngle(oldAngle);
                   // cellMatrix[im][i]->cellType = cellMatrix[im][i+1]->cellType;

                    //double oldAngle = (cellMatrix[im][i-1])->cpCell->angle;
                    //cellMatrix[im][i]->cpCell->angle = cellMatrix[im][i+1]->cpCell->angle;
                    //cellMatrix[im][i]->cellType = cellMatrix[im][i+1]->cellType;
               // }
                cellMatrix[im][in-1]->cpCell->angle = cellMatrix[im][in]->cpCell->angle;
                cellMatrix[im][in-1]->cellType = cellMatrix[im][in]->cellType;
                
            }
            else if(1== bAngle)
            {
                //for (int i = 0; i<im; i++)
                //{
                   // double oldAngle = (cellMatrix[i+1][in])->cpCell->angle;
                   // cellMatrix[i][in]->cpCell->setAngle(oldAngle);
                    //cellMatrix[i][in]->cellType = cellMatrix[i+1][in]->cellType;

                    //double oldAngle = (cellMatrix[i-1][in])->cpCell->angle;
                    //cellMatrix[i][in]->cpCell->angle = cellMatrix[i+1][in]->cpCell->angle;
                    //cellMatrix[i][in]->cellType = cellMatrix[i+1][in]->cellType;
                //}
                cellMatrix[im-1][in]->cpCell->angle = cellMatrix[im][in]->cpCell->angle;
                cellMatrix[im-1][in]->cellType = cellMatrix[im][in]->cellType;

                
            }
            else
            {
                errorCount++;
            }
            newAngle = (bAngle == 0) ?  0.0 : (3.1415/2.0);
            newType = (bType == 0) ? bECOLI_MORAN_A:bECOLI_MORAN_B;
        }
        else if(3 == which)
        {
            newAngle = (bAngle == 0) ?  0.0 : (3.1415/2.0);
            newType = (bType == 0) ? bECOLI_MORAN_A:bECOLI_MORAN_B;
        }




        //output to cell class:
        pCell->cellType = newType;
        pCell->cpCell->setAngle(newAngle);


//        dummyCount = (++dummyCount)%(bM*bN);
//        for(int i = 0;i<bM;i++)
//        {
//            for(int j = 0;j<bN;j++)
//            {
////                //dummy increment
////                if(i*bN +j != dummyCount) continue;

//                jEColi *pCell = cellMatrix[i][j];

//                jEColi_t oldType = pCell->cellType;
//                double oldAngle = pCell->cpCell->angle;

//                // update each cell with logic:
//                double newAngle = ((rand()%2) == 0) ? 0.0:(3.1415/2.0);
//                jEColi_t newType = ((rand()%2) == 0) ? bECOLI_MORAN_A:bECOLI_MORAN_B;

//                //output to cell class:
//                pCell->cellType = newType;
//                cpBodySetAngle(pCell->cpCell->bodyA, newAngle);
//                cpBodySetAngle(pCell->cpCell->bodyB, newAngle);
//            }
//        }

//                // update each cell
//        std::list<jEColi *>::iterator jj;
//        int counter = 0;
//        for ( jj=jcells.begin(); jj!=jcells.end(); jj++ )
//        {
//            jEColi *pCell = (*jj);
//            if (NULL == pCell) continue;
//            counter++;

//            //17Mar.2016:  add a "stuck cell"
//            if(ECOLITYPE_STUCKCELL == (*jj)->cellType) continue;

//            int status      = pCell->jupdate();
//            // check for divisions
//            jEColi * d = (*jj)->jdivide();
//            if ( d != NULL )
//            {
//                jadd_cell ( d );
//                int id = (*jj)->get_id();
//                if(id < SPEED_SIZE)
//                {
//                    dividedCell[id] = 1;
//                }
//            }
//        }//end for

        //jMODS:
//        habitat->stepSimulation(get_chip_dt());
        t += get_sim_dt();
        if ( print_rate > 0 && step % print_rate == 0 )
            print();
        step++;
    }
    else
    {
        emit_message ( "Population limit reached. Increase the population limit via the Simulation menu, or by setting the parameter \"population_max\" in your gro program." );
        set_stop_flag(true);
    }
}
void World::add_signal ( Signal * s ) {

    signal_list.push_back ( s );

}
float World::get_signal_value ( Cell * c, int i ) {

    float s = c->get_size() / 3.0, a = c->get_theta();

    return (
                signal_list[i]->get ( (float) ( c->get_x() ), (float) ( c->get_y() ) ) +
                signal_list[i]->get ( (float) ( c->get_x() + s * cos ( a ) ), (float) ( c->get_y() + s * sin ( a ) ) ) +
                signal_list[i]->get ( (float) ( c->get_x() - s * cos ( a ) ), (float) ( c->get_y() - s * sin ( a ) ) )
                ) / 3.0;

}
void World::emit_signal ( Cell * c, int i, float ds ) {

    signal_list[i]->inc (  c->get_x(), c->get_y(), ds );

}
void World::absorb_signal ( Cell * c, int i, float ds ) {

    signal_list[i]->dec (  c->get_x(), c->get_y(), ds );

}
std::vector< std::vector<float> > * World::get_signal_matrix ( int i ) {

  return signal_list[i]->get_signal_matrix();

}
void World::print ( void ) {

    std::list<Cell *>::iterator j;
    int i;

    for ( j=population->begin(); j!=population->end(); j++ ) {

        printf ( "%d, %d, ", (*j)->get_id(), step );

        for ( i=0; i<MAX_REP_NUM-1; i++ )
            printf ( "%f, ", (*j)->get_fluorescence ( i ) );

        printf ( "%f\n", (*j)->get_fluorescence ( MAX_REP_NUM-1 ) );

    }

}
bool World::snapshot ( const char * path ) {

#ifndef NOGUI
    return calling_thread->snapshot(path);
#endif

}
#define NUM_BINS 12
static int bins[NUM_BINS+1];
static char histbuf[100];
int asd = 0;
void World::histogram ( float x, float y, float width, float height, int channel ) {

    int max_freq;
    std::list<Cell *>::iterator i;
    int j;
    float val;

    for ( j=0; j<NUM_BINS; j++ )
        bins[j]=0;

    max_val = 0.01;
    for ( i=population->begin(); i!=population->end(); i++ ) {
        val = (*i)->get_rep ( channel ) / (*i)->get_size();
        if ( val > max_val )
            max_val = val;
    }

    for ( i=population->begin(); i!=population->end(); i++ ) {
        val = (*i)->get_rep ( channel ) / (*i)->get_size();
        bins[ (int) ( ( NUM_BINS * val) / max_val )]++;
    }

    max_freq = 1;
    for ( j=0; j<NUM_BINS; j++ ) {
        if ( bins[j] > max_freq )
            max_freq = bins[j];
    }
#if PORTED_TO_QT
    if ( channel == 0 )
        glColor3f ( 0.5f, 0.8f, 0.5f);
    else if ( channel == 1 )
        glColor3f ( 0.8f, 0.5f, 0.5f);
    else if ( channel == 2 )
        glColor3f ( 0.8f, 0.8f, 0.5f);

    glBegin(GL_QUADS);
    for ( j=0; j<NUM_BINS; j++ ) {
        glVertex2f ( x + width * (j+0.1) / NUM_BINS, y  );
        glVertex2f ( x + width * (j+0.1) / NUM_BINS, y + height * bins[j] / max_freq );
        glVertex2f ( x + width * (j+1) / NUM_BINS, y + height * bins[j] / max_freq );
        glVertex2f ( x + width * (j+1) / NUM_BINS, y );
    }
    glEnd();

    glLineWidth(1.0f);

    glBegin(GL_LINE_STRIP);
    glVertex2f ( x, y+height );
    glVertex2f ( x, y );
    glVertex2f ( x+width, y );
    glEnd();

    glBegin(GL_LINES);
    glVertex2f ( x+width, y ); glVertex2f ( x+width, y-4 );
    glVertex2f ( x, y+height ); glVertex2f ( x-4, y+height );
    glEnd();

    sprintf ( histbuf, "Channel %d", channel );
    drawString ( x, y + height + 8, histbuf );

    sprintf ( histbuf, "%d", max_freq );
    drawString ( x-24, y + height-4, histbuf );

    sprintf ( histbuf, "%.2f", max_val );
    drawString ( x+width-8, y - 16, histbuf );
#endif

}
void World::scatter ( float x, float y, float width, float height, int channel1, int channel2 ) {

    float max_val = 0.0;
    std::list<Cell *>::iterator i;

    for ( i=population->begin(); i!=population->end(); i++ ) {
        if ( (*i)->get_rep ( channel1 ) / (*i)->get_size() > max_val )
            max_val = (*i)->get_rep ( channel1 ) / (*i)->get_size();
        if ( (*i)->get_rep ( channel2 ) / (*i)->get_size() > max_val )
            max_val = (*i)->get_rep ( channel2 ) / (*i)->get_size();
    }
#if PORTED_TO_QT
    glLineWidth(1.0f);
    glColor3f ( 0.8f, 0.8f, 0.8f);

    glBegin(GL_LINE_STRIP);
    glVertex2f ( x, y+height );
    glVertex2f ( x, y );
    glVertex2f ( x+width, y );
    glEnd();

    glBegin(GL_LINES);
    glVertex2f ( x+width, y ); glVertex2f ( x+width, y-4 );
    glVertex2f ( x, y+height ); glVertex2f ( x-4, y+height );
    glEnd();

    glPointSize ( 2.0 );

    glBegin(GL_POINTS);
    for ( i=population->begin(); i!=population->end(); i++ ) {
        glVertex2f ( x + ( (*i)->get_rep ( channel1 ) / (*i)->get_size() ) * width / max_val,
                     y + ( (*i)->get_rep ( channel2 ) / (*i)->get_size() ) * height / max_val );
    }
    glEnd();

    sprintf ( histbuf, "Channel %d vs. %d", channel1, channel2 );
    drawString ( x, y + height + 8, histbuf );
#endif
}
void World::select_cells ( int x1, int y1, int x2, int y2 ) {

    std::list<Cell *>::iterator i;

    int
            X1 = (1/zoom)*min ( x1, x2 ),
            X2 = (1/zoom)*max ( x1, x2 ),
            Y1 = (1/zoom)*min ( y1, y2 ),
            Y2 = (1/zoom)*max ( y1, y2 );

    for ( i=population->begin(); i!=population->end(); i++ ) {
        if ( X1-5 <= (*i)->get_x() && (*i)->get_x() <= X2+5 && Y1-5 <= (*i)->get_y() && (*i)->get_y() <= Y2+10 ) {
            (*i)->select();
        }
    }

}
void World::deselect_all_cells ( void ) {

    std::list<Cell *>::iterator i;

    for ( i=population->begin(); i!=population->end(); i++ ) {
        (*i)->deselect();
    }

}
void World::dump ( FILE * fp ) {

    std::list<Cell *>::iterator i;

    fprintf ( fp, "id, x, y, theta, volume, gfp, rfp, yfp, cfp\n" );

    for ( i=population->begin(); i!=population->end(); i++ ) {
        fprintf ( fp, "%d, %f, %f, %f, %f, %d, %d, %d, %d\n",
                  (*i)->get_id(), (*i)->get_x(), (*i)->get_y(), (*i)->get_theta(), (*i)->get_volume(),
                  (*i)->get_rep(GFP), (*i)->get_rep(RFP), (*i)->get_rep(YFP), (*i)->get_rep(CFP) );
    }

}
void World::add_barrier ( float x1, float y1, float x2, float y2 ) {

    Barrier * b = new Barrier;

    b->x1 = x1;
    b->y1 = y1;
    b->x2 = x2;
    b->y2 = y2;

    barriers->push_back( *b );

}
