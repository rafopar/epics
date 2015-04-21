#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <iostream>
#include <TGraph.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TApplication.h>

using namespace std;

TH1D *Graph2Hist(TGraph *, double scale = 1.);
double* Calc_abalpha(double, double, double);

bool Fit_2H02A(TGraph *, string);

int main( int argc, char **argv)
{
  TApplication *app1 = new TApplication("", 0, NULL);
  
  if( argc ==1 )
    {
      cout<<"PLease specify at least the harp_name"<<endl;
      cout<<"The program is exiting"<<endl;
      exit(1);
    }
  
  const int n_counters = 15;
  string fname = argv[1];
  cout<<"File name is "<<fname<<endl;
  
  TCanvas *c1 = new TCanvas("c1", "", 900, 900);

  //TGraph *gr1 = new TGraph( "../harp_tagger/harp_tagger_04-20-15_18:26:21.txt", "%lg %*s %*s %*s %*s %*s %lg");
  TGraph *gr1 = new TGraph( "harp_2H02A_04-19-15_09:21:56.txt", "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  //TGraph *gr1 = new TGraph( "harp_2H02A_04-19-15_12:36:21.txt", "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");

  string all_harp_dir = "/home/epics/DATA/HARP_SCANS";
 
  string counter_names_[n_counters] = {"FCup", "Upstream_Left", "Upstream Right", "tagger_left", "tagger_right", "tagger_top", 
				       "downstream_left", "downstream_right",  "downstream_top", "downstream_bottom", 
				       "HPS-Left", "HPS-Right", "HPS-T", "HPS_SC", "empty"};

  TGraph *gr_[n_counters];
  
  gr_[0] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %lg");
  gr_[1] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %lg");
  gr_[2] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %lg");
  gr_[3] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %lg");
  gr_[4] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %lg");
  gr_[5] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %lg");
  gr_[6] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[7] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[8] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[9] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[10] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[11] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[12] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[13] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");
  gr_[14] = new TGraph(Form("%s/%s", all_harp_dir.c_str(), fname.c_str() ), "%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg");

  c1->Divide(1, 14);
  for( int i = 0; i < n_counters - 1; i++ )
    {
      c1->cd(i+1)->SetTopMargin(0.);
      gr_[i]->SetMarkerStyle(23);
      gr_[i]->SetMarkerSize(0.2);
      gr_[i]->SetMarkerColor(2);
      gr_[i]->Draw("AP");
    }

//   gr1->SetMarkerStyle(23);
//   gr1->SetMarkerColor(2);
//   gr1->Draw("AP");

//   TH1D *h_gr1 = (TH1D*)Graph2Hist(gr1);
//   h_gr1->Draw("Same");
  
//   TSpectrum *sp1 = new TSpectrum();
  
//   sp1->Search(h_gr1, 5., "", 0.2);
//   h_gr1->Draw();

//   h_gr1->Sumw2();
  
//   float *peak_val = sp1->GetPositionY();
//   float *pos = sp1->GetPositionX();
//   int n_peaks = sp1->GetNPeaks();
  
//   int n_for_average = 10;
//   double bgr_average;
//   for( int i = 0; i < n_for_average; i++ )
//     {
//       bgr_average = bgr_average + h_gr1->GetBinContent(i+5);
//     }
//   bgr_average = bgr_average/double(n_for_average);
  
//   TF1 *f_GPol3_[n_peaks];
  
//   if( n_peaks == 3 )
//     {
//       for( int i = 0; i < n_peaks; i++ )
// 	{
// 	  f_GPol3_[i] = new TF1(Form("f_GPol3_%d", i), "[0]*TMath::Gaus(x, [1], [2]) + [3]");
// 	  double mean_x = pos[i];
// 	  f_GPol3_[i]->SetRange(mean_x - 0.5, mean_x + 0.5);
// 	  f_GPol3_[i]->SetParameters(peak_val[i] - bgr_average, pos[i], 0.01, bgr_average);
// 	  h_gr1->Fit(f_GPol3_[i], "+MeV", "", mean_x - 0.25, mean_x + 0.25);
// 	  cout<<"Bgr Average is "<<bgr_average<<endl;
// 	  //f_GPol3_[i]->DrawCopy("Same");
// 	}
//     }
  
  Fit_2H02A(gr_[10], counter_names_[10]);
  Fit_2H02A(gr_[11], counter_names_[11]);
  Fit_2H02A(gr_[12], counter_names_[12]);
  Fit_2H02A(gr_[13], counter_names_[13]);
  
  app1->Run();
}

TH1D *Graph2Hist(TGraph * gr, double scale)
{
  gr->Sort();
  int n_points = gr->GetN();
  double *x_axis_tmp = gr->GetX();
  
  cout<<"n_points before "<<n_points<<endl;
  for( int i = 1; i < n_points; i++ )
    {
      cout<<"     ================= Hop ============== x is decreasing ====== i = "<<i<<"    x = "<<x_axis_tmp[i]<<endl;
      if( x_axis_tmp[i] <= x_axis_tmp[i - 1] )
	{
	  gr->RemovePoint(i);
	}
    }

  n_points = gr->GetN();
  double *x_axis = gr->GetX();
  double *y_axis = gr->GetY();

  for( int i = 0; i < n_points; i++ )
    {
      x_axis[i] = scale*x_axis[i];
    }

  double x_varbins_[n_points + 1];
  double y_varbins_[n_points + 1];
  
  cout<<"n_points after rearrangement "<<n_points<<endl;
  for( int i = 1; i < n_points; i++ )
    {
      x_varbins_[i] = (x_axis[i-1] + x_axis[i])/2.;
      //cout<<"x_varbins ["<<i<<"] = "<<x_varbins_[i]<<endl;
      
      if( x_varbins_[i] <= x_varbins_[i - 1] )
	{
	  cout<<"  ================= Hop ============== x is decreasing ====== i = "<<i<<"    x = "<<x_varbins_[i]<<endl;
	}
    }

  x_varbins_[0] = (x_axis[0] - (x_axis[1] - x_axis[0])/2.);
  x_varbins_[n_points] = (x_axis[n_points - 1] + (x_axis[n_points - 1] - x_axis[n_points - 2])/2.);
  
  cout<<"x_varbins ["<<0<<"] = "<<x_varbins_[0]<<endl;
  cout<<"x_varbins ["<<n_points<<"] = "<<x_varbins_[n_points]<<endl;

  double x_low, y_low;
  double x_high, y_high;
  double bin_width = (x_high - x_low)/double(n_points - 1);

  gr->GetPoint(0, x_low, y_low);
  gr->GetPoint(n_points - 1, x_high, y_high);
  TH1D *h_gr = new TH1D("h_gr", "", n_points , x_varbins_ );

  //cout<<"xlow is "<<x_low<<"    x_high is "<<x_high<<endl;

  h_gr->FillN(n_points, x_axis, y_axis);
  h_gr->SetStats(0);
  cout<<"NbinsX = "<<h_gr->GetNbinsX()<<endl;
  
  return h_gr;
}


bool Fit_2H02A(TGraph *gr, string counter_name)
{
  TLatex *lat1 = new TLatex();
  lat1->SetNDC();
  TCanvas *c1 = new TCanvas(Form("%s", counter_name.c_str()), "", 600, 1100);
  TH1D *h_gr = (TH1D*)Graph2Hist(gr, 10.); // 10 is because 2H02A has motor position im cm, it should be converted into mm
  h_gr->SetTitle("; motor pos (mm)");
  h_gr->Sumw2();

  TSpectrum *sp1 = new TSpectrum();
  
  sp1->Search(h_gr, 15., "", 0.2);

  float *peak_val_tmp = sp1->GetPositionY();
  float *pos_tmp = sp1->GetPositionX();
  const int n_peaks = sp1->GetNPeaks();
  
  int n_for_average = 10;
  double bgr_average;
  for( int i = 0; i < n_for_average; i++ )
    {
      bgr_average = bgr_average + h_gr->GetBinContent(i + 5);
    }
  bgr_average = bgr_average/double(n_for_average);
  
  TF1 *f_GPol0_[n_peaks];

  cout<<"Npeaks = "<<n_peaks<<endl;

  if( n_peaks == 3 )
    {
      c1->Clear();
      c1->Divide(1, n_peaks);
      TH1D *h_gr_tmp_[n_peaks];
      
      TGraph *gr_peaks = new TGraph(n_peaks, pos_tmp, peak_val_tmp);
      gr_peaks->Sort();
      double *pos = gr_peaks->GetX();
      double *peak_val = gr_peaks->GetY();

      string wire_names_[3] = {"X", "Y", "45 deg"};
      double mean_[3];
      double sigm_[3];
      double bgr_[3];
      double peak_val_[3];

      for( int i = 0; i < n_peaks; i++ )
	{
	  c1->cd(i+1)->SetLogy();
	  h_gr_tmp_[i] = (TH1D*)h_gr->Clone(Form("h_gr_%d", i));
	  f_GPol0_[i] = new TF1(Form("f_GPol0_%d", i), "[0]*TMath::Gaus(x, [1], [2]) + [3]");
	  double mean_x = pos[i];
	  f_GPol0_[i]->SetRange(mean_x - 5., mean_x + 5.);
	  f_GPol0_[i]->SetParameters(peak_val[i] - bgr_average, pos[i], 0.1, bgr_average);
	  h_gr_tmp_[i]->SetAxisRange(mean_x - 5., mean_x + 5.);
	  //h_gr_tmp_[i]->Draw();
	  h_gr_tmp_[i]->Fit(f_GPol0_[i], "+MeV", "", mean_x - 4.5, mean_x + 4.5);
	  if( i < 2 ) // wires X and Y
	    {
	      mean_[i] = f_GPol0_[i]->GetParameter(1)/sqrt(2.);
	      sigm_[i] = f_GPol0_[i]->GetParameter(2)/sqrt(2.);
	    }
	  else if(i == 2)
	    {
	      mean_[i] = f_GPol0_[i]->GetParameter(1);
	      sigm_[i] = f_GPol0_[i]->GetParameter(2);
	    }
	  bgr_[i] = f_GPol0_[i]->GetParameter(3);
	  peak_val_[i] = f_GPol0_[i]->GetParameter(0);
	  
	  lat1->DrawLatex(0.25, 0.91, Form("Harp: 2H02A   Counter: %s  Wire %s ", counter_name.c_str(), wire_names_[i].c_str()));
	  lat1->DrawLatex(0.12, 0.85, Form("#mu = %1.4f mm", mean_[i]));
	  lat1->DrawLatex(0.12, 0.80, Form("#sigma = %1.4f mm", sigm_[i]));
	  lat1->DrawLatex(0.12, 0.75, Form("peak_val = %1.0f", peak_val[i]));
	  lat1->DrawLatex(0.12, 0.70, Form("bgr_prak = %1.1e", bgr_[i]/peak_val[i]));
	  cout<<"Bgr Average is "<<bgr_average<<endl;
	  //f_GPol0_[i]->DrawCopy("Same");
	}

      double *alpha_a_b = Calc_abalpha(sigm_[2], sigm_[0], sigm_[1]);
      
      double alpha = alpha_a_b[0];
      double aa = alpha_a_b[1];
      double bb = alpha_a_b[2];

      c1->cd(1);
      lat1->DrawLatex(0.7, 0.85, Form("#alpha = %1.2f deg", alpha));
      lat1->DrawLatex(0.7, 0.8, Form("a = %1.2f", aa));
      lat1->DrawLatex(0.7, 0.75, Form("b = %1.2f", bb));

      return true;
    }
  else
    {
      h_gr->Draw();
      lat1->DrawLatex(0.15, 0.91, Form("Harp: 2H02A    counter %s", counter_name.c_str()));
      return false;
    }

}


double* Calc_abalpha(double sigm45, double sigm_x, double sigm_y)
{
  const double radian = 57.2957795130823229;
  const double PI = 3.14159265358979312;
  double A = sigm_x*sigm_x;
  double B = sigm_y*sigm_y;
  double C = sigm45*sigm45;

  double alpha = 0.5*atan((2*C-A-B)/(B-A));
  if( A < B )
    {
      alpha = alpha + PI;
    }
  double a;
  double b;
  
  if( 0.5*(A+B+(A-B)/cos(2*alpha))<0 )
    {
      b = -1.;
    }
  else
    {
      b = sqrt(0.5*( A+B-(A-B)/cos(2*alpha)));
    }

  if( 0.5*(A+B+(A-B)/cos(2*alpha))<0 )
    {
      a = -1.;
    }
  else
    {
      a = sqrt(0.5*(A+B+(A-B)/cos(2*alpha)));
    }
  
  alpha = alpha*radian;
  double alpha_a_b[3];
  alpha_a_b[0] = alpha;
  alpha_a_b[1] = a;
  alpha_a_b[2] = b;
  double *ret_values = alpha_a_b;
  
  return ret_values;
 }
