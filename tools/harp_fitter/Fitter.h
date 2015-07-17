#include <TQObject.h> 
#include <RQ_OBJECT.h>
#include <TGFileDialog.h>

class TF1;
class TH1D;
class TGWindow; 
class TGMainFrame; 
class TRootEmbeddedCanvas;
class TGTextEntry;
class TGTextEdit;
class TGComboBox;
class TGraph;
class TGNumberEntry;

class Fitter { 
   RQ_OBJECT("Fitter")
private: 
  TGMainFrame         *fMain;
  TRootEmbeddedCanvas *fEcanvas;
  TGTextEntry *par0;
  TGFileInfo file_info;
  static const int n_counters = 15;
  static const std::string all_harps_dir;
  void InitData( std::string );
  TGraph *gr_[n_counters];
  std::string counter_names_[n_counters];
  std::string file_name;
  std::string harp_name;
  TGComboBox *counters_box;
  TGTextEdit *comments;
  TGMainFrame *fMain_log;
  bool fit_2c21;
  bool fit_tagger;
  bool fit_2H02A;
  bool preview_mode;
  TGNumberEntry *First_peak_bgr, *First_peak_bgr_min, *First_peak_bgr_max;
  TGNumberEntry *First_peak_A, *First_peak_A_min, *First_peak_A_max;
  TGNumberEntry *First_peak_mean, *First_peak_mean_min, *First_peak_mean_max;
  TGNumberEntry *First_peak_sigm, *First_peak_sigm_min, *First_peak_sigm_max;
  TGNumberEntry *First_peak_range_min, *First_peak_range_max;
  TGNumberEntry *Second_peak_bgr, *Second_peak_bgr_min, *Second_peak_bgr_max;
  TGNumberEntry *Second_peak_A, *Second_peak_A_min, *Second_peak_A_max;
  TGNumberEntry *Second_peak_mean, *Second_peak_mean_min, *Second_peak_mean_max;
  TGNumberEntry *Second_peak_sigm, *Second_peak_sigm_min, *Second_peak_sigm_max;
  TGNumberEntry *Second_peak_range_min, *Second_peak_range_max;
  TGNumberEntry *Third_peak_bgr, *Third_peak_bgr_min, *Third_peak_bgr_max;
  TGNumberEntry *Third_peak_A, *Third_peak_A_min, *Third_peak_A_max;
  TGNumberEntry *Third_peak_mean, *Third_peak_mean_min, *Third_peak_mean_max;
  TGNumberEntry *Third_peak_sigm, *Third_peak_sigm_min, *Third_peak_sigm_max;
  TGNumberEntry *Third_peak_range_min, *Third_peak_range_max;
  
  TF1 *f_1st_peak, *f_2nd_peak, *f_3rd_peak;
  TH1D *h_1st_peak, *h_2nd_peak, *h_3rd_peak;

  double *pars_bgr_1st_peak, *pars_A_1st_peak, *pars_mean_1st_peak, *pars_sigm_1st_peak, *range_1st_peak;
  double *pars_bgr_2nd_peak, *pars_A_2nd_peak, *pars_mean_2nd_peak, *pars_sigm_2nd_peak, *range_2nd_peak;
  double *pars_bgr_3rd_peak, *pars_A_3rd_peak, *pars_mean_3rd_peak, *pars_sigm_3rd_peak, *range_3rd_peak;

public:
  Fitter(const TGWindow *p,UInt_t w,UInt_t h, std::string );
  virtual ~Fitter();
  //void DoDraw();
  void OpenFile();
  void Draw_All_Counters();
  void FitData( bool, bool );
  void GetComments();
  void SubmitToLogbook();
  bool Fit_2c21(TGraph *, std::string );
  bool Search_2c21_peaks(TGraph *);
  void Set_Fit_Pars();
};
