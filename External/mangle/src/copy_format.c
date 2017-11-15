/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <string.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Copy format structure from fmt1 to fmt2, except fmt->nweights (this is used
  as a sort of flag in rdmask as to whether rd_hpix has been called or not).
*/
void copy_format(format *fmt1, format *fmt2)
{
    fmt2->in = fmt1->in;
    fmt2->out = fmt1->out;
    fmt2->skip = fmt1->skip;
    fmt2->end = fmt1->end;
    fmt2->single = fmt1->single;
    fmt2->n = fmt1->n;
    fmt2->nn = fmt1->nn;
    fmt2->innve = fmt1->innve;
    fmt2->outper = fmt1->outper;
    fmt2->outnve = fmt1->outnve;
    fmt2->id = fmt1->id;
    fmt2->newid = fmt1->newid;
    fmt2->pixel = fmt1->pixel;
    fmt2->weight = fmt1->weight;
    fmt2->inunitp = fmt1->inunitp;
    fmt2->outunitp = fmt1->outunitp;
    fmt2->inframe = fmt1->inframe;
    fmt2->outframe = fmt1->outframe;
    fmt2->inunit = fmt1->inunit;
    fmt2->outunit = fmt1->outunit;
    fmt2->outprecision = fmt1->outprecision;
    fmt2->outphase = fmt1->outphase;
    fmt2->azn = fmt1->azn;
    fmt2->eln = fmt1->eln;
    fmt2->azp = fmt1->azp;
    fmt2->trunit = fmt1->trunit;
    fmt2->nweights = fmt2->nweights;
}
