#ifndef _VSMmdbSS_HPP
#define _VSMmdbSS_HPP

#include "hVSMmdbapp.hpp"
#include "hVSMmdbcmd.hpp"
#include "VastSrchUti.hpp"
#include <corelib/ncbiexec.hpp>
#include <corelib/ncbifile.hpp>
#include "mmdbuti.hpp"


void SendSummaryPageText(CCgiContext& ctx, bool VastLink, bool JobDone=true);

#endif
