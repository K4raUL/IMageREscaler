#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <math.h>
#include <random>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    setWindowTitle("IMage REscaler");
    setWindowFlags(Qt::Dialog | Qt::MSWindowsFixedSizeDialogHint);

    start_text = "1. Load Image to display here. Supported extensions:\n"
                 "    BMP, GIF, JPG, JPEG, PNG, PBM, PGM, PPM, XBM, XPM\n"
                 "2. Pick new resolution and RESCALE image\n"
                 "3. RECOLOR image (2 modes):\n"
                 "    3.1. USER mode to manually pick colors from defined preset\n"
                 "    3.2. AUTO mode to generate N colors, corresponding to image palette\n"
                 "4. Save modified Image";
    ui->img_lbl->setText(start_text);

    colorcnt = 5;

    QLineEdit *resle = ui->resw_spin->findChild<QLineEdit*>();
    resle->setReadOnly(0);
    resle->setFocusPolicy(Qt::NoFocus);

    ui->step2->setVisible(0);
    ui->step3->setVisible(0);
    displayDefaultColors();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::clearColors()
{
    delete [] colstr;
    for (unsigned int i = 0; i != COLOR_CNT; i++) delete RBcolor[i];
    delete [] RBcolor;

    clearLayout(ui->color_layout);
}

void MainWindow::clearLayout(QLayout *lo)
{
    QLayoutItem *item;
    while((item = lo->takeAt(0))) {
        delete item->widget();
        delete item;
    }
}

void MainWindow::displayDefaultColors()
{
    clearLayout(ui->color_layout);

    QString pixCol[COLOR_CNT]{"#FFFFFF", "#000000", "#7F7F7F", "#A9A9A9", "#D3D3D3", "#FF0000", "#008000", "#0000FF", "#00FFFF",
                              "#FF00FF", "#FFFF00", "#FFC0CB", "#FFA500", "#90EE90", "#F5F5DC", "#A52A2A", "#ADD8E6", "#B5651D",
                              "#800080", "#87CEEB", "#00FF00", "#FFD700", "#7FFFD4", "#FFE4C4", "#FF7F50", "#FF4500", "#FA8072",
                              "#00008B"};

    colstr = new PIXCOL[COLOR_CNT] {0xFFFFFF , 0x000000 , 0x7F7F7F , 0xA9A9A9 , 0xD3D3D3 , 0xFF0000 , 0x008000 , 0x0000FF , 0x00FFFF ,
                                    0xFF00FF , 0xFFFF00 , 0xFFC0CB , 0xFFA500 , 0x90EE90 , 0xF5F5DC , 0xA52A2A , 0xADD8E6 , 0xB5651D ,
                                    0x800080 , 0x87CEEB , 0x00FF00 , 0xFFD700 , 0x7FFFD4 , 0xFFE4C4 , 0xFF7F50 , 0xFF4500 , 0xFA8072 ,
                                    0x00008B};

    RBcolor = new QCheckBox*[COLOR_CNT];

    for (unsigned int i = 0; i != COLOR_CNT; i++) {
        RBcolor[i] = new QCheckBox("0");
        RBcolor[i]->setStyleSheet("font-weight:bold; background-color:" + pixCol[i] + "; color:" + colorINTtoQSS(contrastColor(colstr[i].col)));
        if (i < 5) RBcolor[i]->setChecked(1);

        connect(RBcolor[i], SIGNAL(clicked(bool)), this, SLOT(updateColorCnt(bool)));
        ui->color_layout->addWidget(RBcolor[i]);
    }
}

void MainWindow::displayModifiedColors(QVector<PIXCOL> clrs, unsigned int *ccnt)
{
    clearLayout(ui->color_layout);

    unsigned int autocolor_cnt = ui->colorcnt_lbl->text().toUInt();
    RBcolor = new QCheckBox*[autocolor_cnt];

    for (unsigned int i = 0; i != autocolor_cnt; i++) {
        RBcolor[i] = new QCheckBox;
        RBcolor[i]->setText(QString::number(ccnt[i]));
        RBcolor[i]->setStyleSheet("background-color:" + colorINTtoQSS(clrs[i].col) + "; color:" + colorINTtoQSS(contrastColor(clrs[i].col)));

        RBcolor[i]->setEnabled(0);
        RBcolor[i]->setChecked(1);

        ui->color_layout->addWidget(RBcolor[i]);
    }
}

void MainWindow::calcClosestGoodRatio()
{
    int srng = 10;      // search range
    int maxGCD = 1, gcd;
    int maxi, maxj;     // final "good" ratio

    for (int i = img_w-srng; i <= img_w; i++) {
        for (int j = img_h-srng; j <= img_h; j++) {
            gcd = GCD(i, j);
            if (gcd > maxGCD) {
                maxi = i;
                maxj = j;
                maxGCD = gcd;
            }
        }
    }

    mod_w = maxi;
    mod_h = maxj;
    maxi /= maxGCD;
    maxj /= maxGCD;

    ui->modres_lbl->setText(QString::number(mod_w) + "x" + QString::number(mod_h));
    ui->modres1_lbl->setText(QString::number(maxi) + "x" + QString::number(maxj));
    ui->resw_spin->setMinimum(maxi);
    ui->resw_spin->setValue(maxi);
    ui->resw_spin->setSingleStep(maxi);
}

void MainWindow::updateColorCnt(bool set)
{
    char delta = set ? 1 : -1;
    colorcnt += delta;
    if (colorcnt)   ui->recolor_btn->setEnabled(1);
    else            ui->recolor_btn->setEnabled(0);
    ui->colorcnt_lbl->setText(QString::number(colorcnt));
}

void MainWindow::on_loadimg_btn_clicked()
{
    QString path = QFileDialog::getOpenFileName(this, tr("Files"), QDir::currentPath(), tr("*.jpg *.png *.bmp *.gif *.jpeg *.pbm *.pgm *.ppm *.xbm *.xpm"));
    if (path.isEmpty()) {
        ui->img_lbl->setText(start_text);
        ui->step2->setVisible(0);
        ui->step3->setVisible(0);
        return;
    }
    filename = path.right(path.size() - path.lastIndexOf('/') - 1);

    qmap.load(path);
    qimg = qmap.toImage();

    img_w = qmap.width();
    img_h = qmap.height();
    imglbl_w = ui->img_lbl->width();
    imglbl_h = ui->img_lbl->height();

    ui->res_lbl->setText(QString::number(img_w) + "x" + QString::number(img_h));
    ui->img_lbl->setPixmap(qmap.scaled(imglbl_w, imglbl_h, Qt::KeepAspectRatio));

    calcClosestGoodRatio();
    ui->step2->setVisible(1);
    ui->step3->setVisible(0);
}


void MainWindow::on_saveimg_btn_clicked()
{
    QString fname = QFileDialog::getSaveFileName(this, tr("Save"), QDir::currentPath(), tr("*.jpg *.png *.bmp *.gif *.jpeg *.pbm *.pgm *.ppm *.xbm *.xpm"));
    QFile f(fname);
    f.open( QIODevice::WriteOnly );
    qimg.save(&f);
    f.close();
}


void MainWindow::on_rescale_btn_clicked()
{
    // get new resolution
    int new_w = ui->resw_spin->value();
    int new_h = ui->resh_lbl->text().toInt();

    qimg = qmap.toImage();
    qimg = qimg.copy((img_w-mod_w)>>1, (img_h-mod_h)>>1, mod_w, mod_h);             // crop original image to get "full adjusted" resolution
    qimg = qimg.scaled(new_w, new_h, Qt::KeepAspectRatio);                          // scale QImage to user resolution (it's size will be reduced)

    qmap_mod = QPixmap::fromImage(qimg);                                            // get new QPixmap from QImage
    ui->img_lbl->setPixmap(qmap_mod.scaled(imglbl_w, imglbl_h, Qt::KeepAspectRatio));   // scale new QPixmap back to ui->img_lbl size and display

    ui->step3->setVisible(1);
}

void MainWindow::on_recolor_btn_clicked()
{
    qimg = qmap_mod.toImage();
    unsigned int *pbit = (unsigned int*)qimg.bits();
    unsigned int pixcnt = qimg.sizeInBytes() >> 2;       // QImage::sizeInBytes() returns 4 times size in pixels ?
    unsigned char colcnt = ui->colorcnt_lbl->text().toInt();

    if (ui->user_btn->isEnabled()) Kmeans(pbit, pixcnt, colcnt);
    else                           Preset(pbit, pixcnt);


    QPixmap qmap1 = qmap.fromImage(qimg);                                            // get new QPixmap from QImage
    ui->img_lbl->setPixmap(qmap1.scaled(imglbl_w, imglbl_h, Qt::KeepAspectRatio));   // scale new QPixmap back to ui->img_lbl size and display
}

// ===================================================== PRESET RECOLOR ====================================================

void MainWindow::Preset(unsigned int *pbit, unsigned int pixcnt)
{
    unsigned int colamount[COLOR_CNT]{0};               // pixels amount of each color

    int i;
    unsigned int minDiff;
    unsigned int diff, imin;

    // get selected colors
    for (unsigned int i = 0; i != COLOR_CNT; i++) {
        if (RBcolor[i]->isChecked()) colstr[i].selected = 1;
        else colstr[i].selected = 0;
    }

    while (pixcnt) {
        pixcnt--;
        minDiff = 0xffffff;

        for (i = 0; i != COLOR_CNT; i++) {
            if (!colstr[i].selected) continue;
            diff = distLABconst((*pbit), colstr[i].L, colstr[i].A, colstr[i].B);
            PIXCOL pb(*pbit);
            //qDebug() << "0x"+QString::number((*pbit)&0xffffff, 16) << " --- " << "0x"+QString::number(colstr[i].col,16) << diff;
            if (diff < 3) {
                imin = i;
                break;
            }
            if (diff < minDiff) {
                minDiff = diff;
                imin = i;
            }
        }
        //qDebug() << "0x"+QString::number((*pbit)&0xffffff, 16) << " --> " << "0x"+QString::number(pixcol[imin],16) << minDiff;
        *pbit = colstr[imin].col;
        colamount[imin]++;
        pbit++;
    }
    for (unsigned int i = 0; i != COLOR_CNT; ++i) RBcolor[i]->setText(QString::number(colamount[i]));
}

// ======================================================== K-MEANS ========================================================

// Defines initial centroids
void MainWindow::Kmeans_init(QVector<PIXCOL> *cntrs, unsigned int *pbit, unsigned int pcnt, unsigned char ccnt) {
    // v1. uniform distribution in RGB-space
    //int step = (1<<24)/(ccnt+1);
    //for (unsigned int i = 1; i != ccnt+1; ++i) cntrs->append(PIXCOL(i*step));

    // v2. uniform distribution in Lab-space
    // ....

    // v3. random pick in RGB-space (! possible duplicates !)
    //std::mt19937 mt(time(NULL));
    //for (unsigned int i = 0; i != ccnt; ++i) cntrs->append(PIXCOL(mt() % (1<<24)));

    // v4. Kmeans++ initial centroids (1st is random)
    std::mt19937 mt(time(NULL));
    double minDiff, maxDiff, diff;
    unsigned int pixs, maxD2;

    cntrs->append(PIXCOL( (*(pbit + mt()%pcnt)) & 0xffffff ));     // first centroid picked randomly from one of pixels
    for (unsigned char i = 0; i != ccnt-1; ++i) {
        pixs = 0;
        maxD2 = 0;
        maxDiff = 0;

        while (pixs != pcnt) {
            minDiff = 0xffffff;

            // Find closest centroid; (i+1) equals cntrs->size()
            for (unsigned char j = 0; j != i+1; ++j) {
                diff = distLABconst(*(pbit+pixs), cntrs->at(j).L, cntrs->at(j).A, cntrs->at(j).B);
                if (diff < minDiff) minDiff = diff;
            }

            if (minDiff > maxDiff) {
                maxDiff = minDiff;
                maxD2 = *(pbit+pixs);
            }

            pixs++;
        }

        cntrs->append(PIXCOL(maxD2 & 0xffffff));
    }
}

double MainWindow::Kmeans_step(QVector<PIXCOL> *cntrs, unsigned int *pbit, unsigned int pcnt, unsigned char ccnt)
{
    double sumDiff[ccnt];                                   // summ distance of each cluster pixel to its centroid
    unsigned int klas_col[ccnt];                            // amount of pixels, assigned to each cluster/centroid
    double newCML[ccnt], newCMA[ccnt], newCMB[ccnt];        // summ LAB components of pixels for each cluster, to calculate average (for new centroids)
    //unsigned int newCMr[ccnt], newCMg[ccnt], newCMb[ccnt];  // summ RGB components of pixels for each cluster, to calculate average (for new centroids)

    memset(sumDiff, 0, sizeof(double)*ccnt);
    memset(klas_col, 0, sizeof(unsigned int)*ccnt);
    memset(newCML, 0, sizeof(double)*ccnt);
    memset(newCMA, 0, sizeof(double)*ccnt);
    memset(newCMB, 0, sizeof(double)*ccnt);
    //memset(newCMr, 0, sizeof(unsigned int)*ccnt);
    //memset(newCMg, 0, sizeof(unsigned int)*ccnt);
    //memset(newCMb, 0, sizeof(unsigned int)*ccnt);


    double minDiff, diff;
    double sumDev = 0;
    unsigned int pixs = 0;
    unsigned char i, mini;
    //unsigned char r, g, b;

    while(pixs != pcnt) {
        minDiff = 0xffffff;
        PIXCOL npix(*(pbit+pixs));

        // Assign each pixel to closest centroid
        for (i = 0; i != ccnt; ++i) {
            diff = distLAB2(npix.L, npix.A, npix.B, cntrs->at(i).L, cntrs->at(i).A, cntrs->at(i).B);        // simple euclidean distance squared
            //diff = distRGB(*(pbit+pixs), cntrs->at(i).col);
            if (diff < minDiff) {
                minDiff = diff;
                mini = i;
            }
        }
        sumDiff[mini] += minDiff;
        klas_col[mini]++;

        // RGB components
        //r = *(pbit+pixs) >> 16; g = *(pbit+pixs) >> 8; b = *(pbit+pixs);
        //newCMr[mini] += r;
        //newCMg[mini] += g;
        //newCMb[mini] += b;

        // LAB components
        newCML[mini] += npix.L;
        newCMA[mini] += npix.A;
        newCMB[mini] += npix.B;

        pixs++;
    }

    for (i = 0; i != ccnt; ++i) {
        if (klas_col[i] == 0) continue;     // no pixels were assigned to this centroid
        sumDiff[i] /= klas_col[i];

        // average in RGB space (! integer division !)
        //newCMr[i] /= klas_col[i];
        //newCMg[i] /= klas_col[i];
        //newCMb[i] /= klas_col[i];

        // average in LAB space
        newCML[i] /= klas_col[i];
        newCMA[i] /= klas_col[i];
        newCMB[i] /= klas_col[i];

        //PIXCOL ncol((newCMr[i]<<16) | (newCMg[i] << 8) | newCMb[i]);
        PIXCOL ncol(LABtoRGB(newCML[i], newCMA[i], newCMB[i]));
        cntrs->replace(i, ncol);
        sumDev += sumDiff[i];
    }

    return sumDev;
}

void MainWindow::Kmeans(unsigned int *pbit, unsigned int pixcnt, unsigned char colcnt)
{
    QVector<PIXCOL> centers;
    Kmeans_init(&centers, pbit, pixcnt, colcnt);

    // Kmeans iterations
    double sumDev = 0, sumDevmin;
    do {
        sumDevmin = sumDev;
        sumDev = Kmeans_step(&centers, pbit, pixcnt, colcnt);
    } while (fabs(sumDevmin - sumDev) > 1.);

    // apply new centers colors (modifying *pbit)
    double diff, minDiff;
    unsigned int klas_col[colcnt];
    int pixs = pixcnt;
    unsigned char i, mini;

    memset(klas_col, 0, sizeof(unsigned int)*colcnt);

    while(pixs) {
        pixs--;
        minDiff = 0xffffff;

        for (i = 0; i != colcnt; ++i) {
            diff = distLABconst((*pbit), centers[i].L, centers[i].A, centers[i].B);
            if (diff < minDiff) {
                minDiff = diff;
                mini = i;
            }
        }
        *pbit = centers[mini].col;
        klas_col[mini]++;
        pbit++;
    }

    // redraw colors
    displayModifiedColors(centers, klas_col);
}

// ======================================================== FUNCTIONS ========================================================

double distLAB(int c1, int c2)
{
    double L1, A1, B1, L2, A2, B2;
    RGBtoLAB(c1, &L1, &A1, &B1);
    RGBtoLAB(c2, &L2, &A2, &B2);

    return (L1-L2)*(L1-L2) + (A1-A2)*(A1-A2) + (B1-B2)*(B1-B2);
}

double distLABconst(int c1, double L, double A, double B)
{
    double L1, A1, B1;
    RGBtoLAB(c1, &L1, &A1, &B1);

    return (L1-L)*(L1-L) + (A1-A)*(A1-A) + (B1-B)*(B1-B);
}

double distLAB2(double L1, double A1, double B1, double L2, double A2, double B2)
{
    return (L1-L2)*(L1-L2) + (A1-A2)*(A1-A2) + (B1-B2)*(B1-B2);
}

double distRGB(int c1, int c2)
{
    unsigned char r1 = c1>>16, g1 = c1>>8, b1 = c1;
    unsigned char r2 = c2>>16, g2 = c2>>8, b2 = c2;

    double rmean = (r1 + r2)/512;                   // in fact rmean/256
    int dr = r1-r2, dg = g1-g2, db = b1-b2;

    return (2 + rmean)*dr*dr + 4*dg*dg + (2.99609375 - rmean)*db*db;        // 2 + 255/256 = 2.99609375;
}

void RGBtoXYZ(int rgb, double *x, double *y, double *z)
{
    unsigned char r = rgb>>16;
    unsigned char g = rgb>>8;
    unsigned char b = rgb;
    double var_R = r/255.;
    double var_G = g/255.;
    double var_B = b/255.;

    if ( var_R > 0.04045 ) var_R = pow((var_R + 0.055)/1.055, 2.4);
    else                   var_R = var_R / 12.92;
    if ( var_G > 0.04045 ) var_G = pow((var_G + 0.055)/1.055, 2.4);
    else                   var_G = var_G / 12.92;
    if ( var_B > 0.04045 ) var_B = pow((var_B + 0.055)/1.055, 2.4);
    else                   var_B = var_B / 12.92;

    *x = var_R * 41.24 + var_G * 35.76 + var_B * 18.05;
    *y = var_R * 21.26 + var_G * 71.52 + var_B * 7.22;
    *z = var_R * 1.93 + var_G * 11.92 + var_B * 95.05;
}

void XYZtoLAB(double x, double y, double z, double *L, double *a, double *b)
{
    double var_X = x / 95.047;
    double var_Y = y / 100.;
    double var_Z = z / 108.883;

    if ( var_X > 0.008856 ) var_X = pow(var_X, 1./3.);
    else                    var_X = ( 7.787 * var_X ) + 0.13793103448;
    if ( var_Y > 0.008856 ) var_Y = pow(var_Y, 1./3.);
    else                    var_Y = ( 7.787 * var_Y ) + 0.13793103448;
    if ( var_Z > 0.008856 ) var_Z = pow(var_Z, 1./3.);
    else                    var_Z = ( 7.787 * var_Z ) + 0.13793103448;

    *L = ( 116 * var_Y ) - 16;
    *a = 500 * ( var_X - var_Y );
    *b = 200 * ( var_Y - var_Z );
}

void RGBtoLAB(int rgb, double *L, double *A, double *B)
{
    unsigned char r = rgb>>16;
    unsigned char g = rgb>>8;
    unsigned char b = rgb;
    double var_R = r/255.;
    double var_G = g/255.;
    double var_B = b/255.;

    if ( var_R > 0.04045 ) var_R = pow((var_R + 0.055)/1.055, 2.4);
    else                   var_R = var_R / 12.92;
    if ( var_G > 0.04045 ) var_G = pow((var_G + 0.055)/1.055, 2.4);
    else                   var_G = var_G / 12.92;
    if ( var_B > 0.04045 ) var_B = pow((var_B + 0.055)/1.055, 2.4);
    else                   var_B = var_B / 12.92;

    double var_X = (var_R*0.4338906015 + var_G*0.376235 + var_B*0.189906046);
    double var_Y = (var_R*0.2126 + var_G*0.7152 + var_B*0.0722);
    double var_Z = (var_R*0.017725448 + var_G*0.109475308 + var_B*0.872955374);

    if ( var_X > 0.008856 ) var_X = pow(var_X, 1./3.);
    else                    var_X = ( 7.787 * var_X ) + 0.13793103448;
    if ( var_Y > 0.008856 ) var_Y = pow(var_Y, 1./3.);
    else                    var_Y = ( 7.787 * var_Y ) + 0.13793103448;
    if ( var_Z > 0.008856 ) var_Z = pow(var_Z, 1./3.);
    else                    var_Z = ( 7.787 * var_Z ) + 0.13793103448;

    *L = 116*var_Y - 16;
    *A = 500*( var_X - var_Y );
    *B = 200*( var_Y - var_Z );
}

void LABtoXYZ(double L, double a, double b, double *x, double *y, double *z)
{
    double var_Y = ( L + 16 ) / 116;
    double var_X = 0.002*a + var_Y;
    double var_Z = var_Y - 0.005*b;

    if ( var_Y  > 0.20689303442 ) var_Y = var_Y*var_Y*var_Y;
    else                       var_Y = (var_Y - 0.13793103448) / 7.787;
    if ( var_X  > 0.20689303442 ) var_X = var_X*var_X*var_X;
    else                       var_X = (var_X - 0.13793103448) / 7.787;
    if ( var_Z  > 0.20689303442 ) var_Z = var_Z*var_Z*var_Z;
    else                       var_Z = (var_Z - 0.13793103448) / 7.787;

    *x = var_X * 95.047;
    *y = var_Y * 100;
    *z = var_Z * 108.883;
}

void XYZtoRGB(double x, double y, double z, int *rgb)
{
    double var_X = x*0.01;
    double var_Y = y*0.01;
    double var_Z = z*0.01;

    double var_R = var_X *  3.2406 - var_Y * 1.5372 - var_Z * 0.4986;
    double var_G = -var_X * 0.9689 + var_Y * 1.8758 + var_Z * 0.0415;
    double var_B = var_X *  0.0557 - var_Y * 0.2040 + var_Z * 1.0570;

    if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R, 1./2.4) - 0.055;
    else                     var_R = 12.92 * var_R;
    if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G, 1./2.4) - 0.055;
    else                     var_G = 12.92 * var_G;
    if ( var_B > 0.0031308 ) var_B = 1.055 * pow(var_B, 1./2.4) - 0.055;
    else                     var_B = 12.92 * var_B;

    unsigned char sR = var_R * 255;
    unsigned char sG = var_G * 255;
    unsigned char sB = var_B * 255;

    *rgb = (sR<<16) | (sG<<8) | sB;
}

int LABtoRGB(double L, double A, double B)
{
    double var_Y = ( L + 16 ) / 116;
    double var_X = 0.002*A + var_Y;
    double var_Z = var_Y - 0.005*B;

    if ( var_Y  > 0.20689303442 ) var_Y = var_Y*var_Y*var_Y;
    else                       var_Y = (var_Y - 0.13793103448) / 7.787;
    if ( var_X  > 0.20689303442 ) var_X = var_X*var_X*var_X;
    else                       var_X = (var_X - 0.13793103448) / 7.787;
    if ( var_Z  > 0.20689303442 ) var_Z = var_Z*var_Z*var_Z;
    else                       var_Z = (var_Z - 0.13793103448) / 7.787;

    double var_R =  var_X * 3.080093083 - var_Y * 1.5372 - var_Z * 0.542890638;
    double var_G = -var_X * 0.920910383 + var_Y * 1.8758 + var_Z * 0.045186445;
    double var_B =  var_X * 0.052941179 - var_Y * 0.2040 + var_Z * 1.150893310;

    if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R, 1./2.4) - 0.055;
    else                     var_R = 12.92 * var_R;
    if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G, 1./2.4) - 0.055;
    else                     var_G = 12.92 * var_G;
    if ( var_B > 0.0031308 ) var_B = 1.055 * pow(var_B, 1./2.4) - 0.055;
    else                     var_B = 12.92 * var_B;

    unsigned char sR = var_R * 255;
    unsigned char sG = var_G * 255;
    unsigned char sB = var_B * 255;

    return (sR<<16) | (sG<<8) | sB;
}

int GCD(int n1, int n2)
{
    while(n1 != n2) {
        if (n1 > n2) n1 -= n2;
        else         n2 -= n1;
    }
    return n1;
}

int invColor(int c)
{
    unsigned char r = c>>16;
    unsigned char g = c>>8;
    unsigned char b = c;
    r = 255-r;
    g = 255-g;
    b = 255-b;

    return (r<<16) | (g<<8) | b;
}

int contrastColor(int c)
{
    unsigned char r = c>>16;
    unsigned char g = c>>8;
    unsigned char b = c;

    double Y = 0.0011765*r + 0.0023137255*g + 0.00043137*b;        // luminance. https://en.wikipedia.org/wiki/HSL_and_HSV#Lightness
    return (Y >= 0.5) ? 0 : (1<<24)-1;
}

QString colorINTtoQSS(int c)
{
    QString qcol = QString::number(c, 16);
    while (qcol.size() != 6) qcol = "0" + qcol;
    return "#" + qcol;
}

unsigned int colorQSStoINT(QString c)
{
    bool ok;
    while (c[0] == '#' || c[0] == '0') c.remove(0, 1);
    return c.toUInt(&ok, 16);
}

// ======================================================== SLOTS ========================================================

void MainWindow::on_resw_spin_valueChanged(int w1)
{
    ui->resh_lbl->setText(QString::number(w1*mod_h/mod_w));
}

void MainWindow::on_colall_btn_clicked()
{
    if (colorcnt < COLOR_CNT){
        colorcnt = COLOR_CNT;
        for (unsigned int i = 0; i != COLOR_CNT; i++) RBcolor[i]->setChecked(1);
        ui->colorcnt_lbl->setText(QString::number(COLOR_CNT));
    }
    else {
        colorcnt = 0;
        for (unsigned int i = 0; i != COLOR_CNT; i++) RBcolor[i]->setChecked(0);
        ui->colorcnt_lbl->setText("0");
        ui->recolor_btn->setEnabled(0);
    }
}

void MainWindow::on_bw_btn_clicked()
{
    for (unsigned int i = 0; i != COLOR_CNT; i++) {
        if (i < 5) RBcolor[i]->setChecked(1);
        else RBcolor[i]->setChecked(0);
    }
    colorcnt = 5;
    ui->colorcnt_lbl->setText("5");
}

void MainWindow::on_auto_btn_clicked()
{
    ui->user_btn->setStyleSheet("");
    ui->auto_btn->setStyleSheet("background-color:lightgreen;");

    ui->auto_btn->setEnabled(0);
    ui->user_btn->setEnabled(1);

    ui->colorcnt_lbl->setEnabled(1);
    ui->colorcnt_lbl->setReadOnly(0);

    ui->label_7->setText("/32");

    ui->colall_btn->setEnabled(0);
    ui->bw_btn->setEnabled(0);

    clearColors();
}

void MainWindow::on_user_btn_clicked()
{
    ui->auto_btn->setStyleSheet("");
    ui->user_btn->setStyleSheet("background-color:lightgreen;");

    ui->user_btn->setEnabled(0);
    ui->auto_btn->setEnabled(1);

    ui->label_7->setText("/28");

    ui->colall_btn->setEnabled(1);
    ui->bw_btn->setEnabled(1);

    ui->colorcnt_lbl->setEnabled(0);
    ui->colorcnt_lbl->setReadOnly(1);
    ui->colorcnt_lbl->setText("5");
    clearLayout(ui->color_layout);
    displayDefaultColors();
}


void MainWindow::on_colorcnt_lbl_editingFinished()
{
    unsigned int ccnt = ui->colorcnt_lbl->text().toUInt();
    if (ccnt > 32) ui->colorcnt_lbl->setText("32");
    else if (ccnt < 2) ui->colorcnt_lbl->setText("2");
}

