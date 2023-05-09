#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
#include <QPixmap>
#include <QLabel>
#include <QCheckBox>

#define COLOR_CNT 28

// ==========================================================================================================================================================================
double distLAB(int c1, int c2);                                 // squared euclidean distance in Lab color space
double distLABconst(int c1, double L, double A, double B);
double distLAB2(double L1, double A1, double B1, double L2, double A2, double B2);
double distRGB(int c1, int c2);                                 // squared euclidean distance approximation in RGB color space

void RGBtoXYZ(int rgb, double *x, double *y, double *z);
void XYZtoLAB(double x, double y, double z, double *L, double *a, double *b);
void RGBtoLAB(int rgb, double *L, double *a, double *b);

void LABtoXYZ(double L, double a, double b, double *x, double *y, double *z);
void XYZtoRGB(double x, double y, double z, int *rgb);
int  LABtoRGB(double L, double A, double B);

int GCD(int n1, int n2);
int invColor(int c);
int contrastColor(int c);
QString colorINTtoQSS(int c);
unsigned int colorQSStoINT(QString c);

typedef struct PIXCOL {
    PIXCOL(int arg) {
        col = arg;
        r = arg>>16;    g = arg>>8;     b = arg;
        RGBtoLAB(arg, &L, &A, &B);
    }
    //PIXCOL(PIXCOL &cpy) {col = cpy.col;}
    double L, A, B;
    int col;                // 24bits (32 with transparency)
    unsigned char r, g, b;
    bool selected = 0;
}PIXCOL;
// ==========================================================================================================================================================================

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void calcClosestGoodRatio();
    void displayDefaultColors();
    void displayModifiedColors(QVector<PIXCOL>, unsigned int *ccnt);

    void clearColors();
    void clearLayout(QLayout *lo);

    void Preset(unsigned int *pbit, unsigned int pixcnt);                           // modifies current qimg.bits() !!
    void Kmeans(unsigned int *pbit, unsigned int pixcnt, unsigned char colcnt);     // modifies current qimg.bits() !!
    void Kmeans_init(QVector<PIXCOL> *cntrs, unsigned int *pbit, unsigned int pcnt, unsigned char ccnt);
    double Kmeans_step(QVector<PIXCOL> *cntrs, unsigned int *pbit, unsigned int pcnt, unsigned char ccnt);

private slots:
    void on_loadimg_btn_clicked();
    void on_saveimg_btn_clicked();

    void on_rescale_btn_clicked();
    void on_recolor_btn_clicked();

    void updateColorCnt(bool);

    void on_resw_spin_valueChanged(int w1);
    void on_colall_btn_clicked();
    void on_bw_btn_clicked();

    void on_auto_btn_clicked();
    void on_user_btn_clicked();

    void on_colorcnt_lbl_editingFinished();

private:
    Ui::MainWindow *ui;
    QImage qimg;
    QPixmap qmap, qmap_mod;

    QCheckBox **RBcolor;
    PIXCOL *colstr;
    unsigned short colorcnt;

    int imglbl_w, imglbl_h;     // image size inside QLabel
    int img_w, img_h;           // real image size
    int mod_w, mod_h;           // adjusted image size

    QString start_text;
    QString filename;


    //bool eventFilter(QObject *object, QEvent *event);
};
#endif // MAINWINDOW_H
