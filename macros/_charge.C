#include "../style/Style.C"
#include "../style/Labels.C"
void _charge()
{

  
  SetQQbarStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(0);

  gStyle->SetTitleX(0.17); 
  gStyle->SetTitleY(0.83);

  // Modify accordingly the inputs in case of using ttbar samples... also the "bbbar" string in the first process line

  gROOT->ProcessLine(".x chargeEfficiency.C (\"/home/irles/WorkArea/BBbar_tests/ntuples/trackrestoring/IDR_vertex_reprocessed/OUTPUT3_eL_l5_valencia2.0_merged.root\",kGray+2, NULL,\"/home/irles/WorkArea/BBbar_tests/ntuples/trackrestoring/IDR_vertex_reprocessed/OUTPUT2_eL_l5_valencia2.0_merged.root\",\"ttbar\")");  
  gROOT->ProcessLine(".x chargeEfficiency.C (\"/home/irles/WorkArea/BBbar_tests/ntuples/trackrestoring/IDR_vertex_reprocessed/OUTPUT5_eL_l5_valencia2.0_merged.root\",kBlue, c1,\"/home/irles/WorkArea/BBbar_tests/ntuples/trackrestoring/IDR_vertex_reprocessed/OUTPUT2_eL_l5_valencia2.0_merged.root\")");


}
