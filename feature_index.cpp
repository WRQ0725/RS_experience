#include "head.h"

//1������RVIָ��
Mat CalRVI(Mat ImgR, Mat ImgNir) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//��ȡͼ��ĸ�
    int nWidth = ImgR.cols;//��ȡͼ��Ŀ�
    int nChannels = ImgR.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImDataNir = (double*)ImgNir.data; //

    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)temp2 / temp1;
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}
//2������NDVIָ��
Mat CalNDVI(Mat ImgR, Mat ImgNir) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//��ȡͼ��ĸ�
    int nWidth = ImgR.cols;//��ȡͼ��Ŀ�
    int nChannels = ImgR.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImDataNir = (double*)ImgNir.data; //
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNDVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            int temp1 = *(pImDataR + i * nWidth + j);
            int temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)(temp2 - temp1) / (double)(temp1 + temp2);
            *(pImDataNDVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

//3������SAVIָ��
Mat CalSAVI(Mat ImgR, Mat ImgNir) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//��ȡͼ��ĸ�
    int nWidth = ImgR.cols;//��ȡͼ��Ŀ�
    int nChannels = ImgR.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImDataNir = (double*)ImgNir.data; //

    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNDVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)(temp2 - temp1) / (double)(temp1 + temp2 + 0.5) * (1 + 0.5);
            *(pImDataNDVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}
//4�������������SAVIָ��
Mat CalSAVI2(Mat ImgR, Mat ImgNir) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//��ȡͼ��ĸ�
    int nWidth = ImgR.cols;//��ȡͼ��Ŀ�
    int nChannels = ImgR.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImDataNir = (double*)ImgNir.data; //
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNDVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (temp2 * 2 + 1 - sqrt(pow((2 * temp2 + 1), 2) - 8 * (temp2 - temp1))) / 2;
            *(pImDataNDVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

//5������TGDVIָ��
Mat CalTGDVI(Mat ImgR, Mat ImgNir,Mat ImgG) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//��ȡͼ��ĸ�
    int nWidth = ImgR.cols;//��ȡͼ��Ŀ�
    int nChannels = ImgR.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImDataNir = (double*)ImgNir.data; 
    double* pImDataG = (double*)ImgG.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNDVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = *(pImDataG + i * nWidth + j);
            double temp4 = (temp2 - temp1)/(0.83-0.66)-(temp1-temp3)/(0.66-0.56);
            *(pImDataNDVI + i * nWidth + j) = temp4;
        }
    }
    return M;
}

//6������NDWIָ��
Mat CalNDWI(Mat ImgG, Mat ImgNir) {
    double* pImDataR = (double*)ImgG.data;
    int nHeight = ImgG.rows;//��ȡͼ��ĸ�
    int nWidth = ImgG.cols;//��ȡͼ��Ŀ�
    int nChannels = ImgG.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImDataNir = (double*)ImgNir.data; //
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNDVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (temp1 - temp2) / (temp1 + temp2);
            *(pImDataNDVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

//7������NWIָ��
Mat CalNWI(Mat Img1, Mat Img4, Mat Img5, Mat Img7) {
    double* pImData[4];
    pImData[0] = (double*)Img1.data;
    int nHeight = Img1.rows;//��ȡͼ��ĸ�
    int nWidth = Img1.cols;//��ȡͼ��Ŀ�
    int nChannels = Img1.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    pImData[1] = (double*)Img4.data;
    pImData[2] = (double*)Img5.data;
    pImData[3] = (double*)Img7.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNWI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp[4];
            for (int k = 0; k < 4; k++) {
                temp[k] = *(pImData[k] + i * nWidth + j);
            }
            double NWI = (temp[0] - temp[1] - temp[2] - temp[3]) / (temp[0] + temp[1] + temp[2] + temp[3]) * 255;
            *(pImDataNWI + i * nWidth + j) = NWI;
        }
    }
    return M;
}

//8��DBI����Ҫ�����ֽ����õغ�ֲ�����ǵ���
Mat CalDBI(Mat Img4, Mat Img7) {
    double* pImData4 = (double*)Img4.data;
    int nHeight = Img4.rows;//��ȡͼ��ĸ�
    int nWidth = Img4.cols;//��ȡͼ��Ŀ�
    int nChannels = Img4.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImData7 = (double*)Img7.data; //
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataDBI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImData4 + i * nWidth + j);
            double temp2 = *(pImData7 + i * nWidth + j);
            double temp3 = temp2 - temp1;
            *(pImDataDBI + i * nWidth + j) = temp3;
        }
    }
    return M;
}
//9��NDBIָ��(��Ҫ��5���λ���7���Σ���������SWIR�����ǲ���λ�������𣬽�����SWIR-2���������Ĳ�����
Mat CalNDBI(Mat Img4, Mat Img7) {
    double* pImData4 = (double*)Img4.data;
    int nHeight = Img4.rows;//��ȡͼ��ĸ�
    int nWidth = Img4.cols;//��ȡͼ��Ŀ�
    int nChannels = Img4.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImData5 = (double*)Img7.data; //
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNDBI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImData4 + i * nWidth + j);
            double temp2 = *(pImData5 + i * nWidth + j);
            double temp3 = (temp2 - temp1) / (temp2 + temp1);
            *(pImDataNDBI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

//10��IBIָ����ͬ�ϣ�����Ϊ����7��
Mat CalIBI(Mat Img2, Mat Img3, Mat Img4, Mat Img7) {
    double* pImData1 = (double*)Img2.data;
    int nHeight = Img2.rows;//��ȡͼ��ĸ�
    int nWidth = Img2.cols;//��ȡͼ��Ŀ�
    int nChannels = Img2.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImData2 = (double*)Img3.data;
    double* pImData3 = (double*)Img4.data;
    double* pImData7 = (double*)Img7.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataIBI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImData1 + i * nWidth + j);
            double temp2 = *(pImData2 + i * nWidth + j);
            double temp3 = *(pImData3 + i * nWidth + j);
            double temp4 = *(pImData7 + i * nWidth + j);
            double IBI = (2 * temp4 / (temp4 + temp3) - (temp3 / (temp3 + temp2)) - (temp1 / (temp1 + temp4))) / (2 * temp4 / (temp4 + temp3) + (temp3 / (temp3 + temp2)) + (temp1 / (temp1 + temp4)));
            *(pImDataIBI + i * nWidth + j) = IBI;
        }
    }
    return M;
}

//�Զ����ͼ�񣨵����Σ�
bool ReadImage(char imageName[], Mat& M) {
    M = imread(imageName, IMREAD_ANYCOLOR);
    if (M.empty())     // �ж��ļ��Ƿ�������     
    {
        fprintf(stderr, "Can not load image %s\n", imageName);
        waitKey(6000);  // �ȴ�6000 ms�󴰿��Զ��ر�   
        return false;
    }
}


//�Ҷ��������죬ʹ������ͼ�нϺõ�Ŀ��Ч��
Mat LineTransform(Mat M, float a, float b) {
    double* pImData = (double*)M.data;
    int nHeight = M.rows;//��ȡͼ��ĸ�
    int nWidth = M.cols;//��ȡͼ��Ŀ�
    int nChannels = M.channels();//��ȡͼ��ͨ����Ŀ BGR
    Mat Mnew = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNew = (double*)Mnew.data;
    for (int k = 0; k < nChannels; k++) {//ͳ�Ƴ�ԭͼ��������С�Ҷȼ�
        double nMingray = 255, nMaxgray = 0;
        for (int i = 0; i < nHeight; i++) {
            for (int j = 0; j < nWidth; j++) {
                double temp = *(pImData + (nWidth * i + j) * nChannels + k);
                if (temp < nMingray) {
                    nMingray = temp;
                }
                if (temp > nMaxgray) {
                    nMaxgray = temp;
                }
            }
        }
        for (int i = 0; i < nHeight; i++) {
            for (int j = 0; j < nWidth; j++) {
                double temp = *(pImData + (nWidth * i + j) * nChannels + k);
                temp = (double)(a + (b - a) / (nMaxgray - nMingray) * (temp - nMingray));

                *(pImDataNew + (nWidth * i + j) * nChannels + k) = temp;
            }
        }
    }
    return Mnew;
}
//ֱ��ͼͳ�Ʒ����ж�ֵ������
bool StateMethod2(Mat& M) {
    unsigned char* pImData = M.data;
    int nHeight = M.rows;//��ȡͼ��ĸ�
    int nWidth = M.cols;//��ȡͼ��Ŀ�
    int nChannels = M.channels();//��ͨ��
    int nThreshold = 0;//��ʼ������ֵ
    int nNewThreshold = 0;//�����µ���ֵ
    int  hist[256];//ͳ�ƻҶȼ�
    memset(hist, 0, sizeof(hist));
    int nTotalState1 = 0, nTotalState2 = 0;
    double dRecState1 = 0, dRecState2 = 0, dGrayAvgState1 = 0, dGrayAvgState2 = 0;
    int nItereationTimes = 0;
    int nMaxGray = 0, nMinGray = 255;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            int temp = *(pImData + i * nWidth + j);
            hist[temp]++;//ͳ��Ƶ��
            if (nMaxGray < temp) {
                nMaxGray = temp;
            }
            if (nMinGray > temp) {
                nMinGray = temp;
            }
        }
    }
    nNewThreshold = (nMaxGray + nMinGray) / 2;
    for (; nThreshold != nNewThreshold && nItereationTimes < 100; nItereationTimes++) {
        nThreshold = nNewThreshold;
        nTotalState1 = 0, nTotalState2 = 0;
        dRecState1 = 0.0, dRecState2 = 0.0;
        for (int k = nMinGray; k < nThreshold; k++) {
            dRecState1 += (double)k * hist[k];
            nTotalState1 += hist[k];
        }
        dGrayAvgState1 = dRecState1 / nTotalState1;
        for (int k = nThreshold; k <= nMaxGray; k++) {
            dRecState2 += (double)k * hist[k];
            nTotalState2 += hist[k];
        }
        dGrayAvgState2 = dRecState2 / nTotalState2;
        nNewThreshold = (int)((dGrayAvgState1 + dGrayAvgState2) / 2);
    }
    //ѭ������ͼ������ֵ���ж�ֵ��
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            if (*(pImData + i * nWidth + j) < nThreshold) {
                *(pImData + i * nWidth + j) = 0;
            }
            else {
                *(pImData + i * nWidth + j) = 255;
            }
        }
    }
    cout << nThreshold;
    return true;
}





void Run_RVI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalRVI(M[2], M[3]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("RVI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_WINTER);
    imwrite("RVI_bin.jpg", Mth);
}

void Run_NDVI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalNDVI(M[2], M[3]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("NDVI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_WINTER);
    imwrite("NDVI_bin.jpg", Mth);
}
void Run_SAVI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalSAVI(M[2], M[3]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("SAVI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_WINTER);
    imwrite("SAVI_bin.jpg", Mth);
}

void Run_SAVI2(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalSAVI2(M[2], M[3]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("SAVI2.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_WINTER);
    imwrite("SAVI2_bin.jpg", Mth);
}

void Run_TGDVI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalTGDVI(M[2], M[3],M[1]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("TGDVI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_WINTER);
    imwrite("TGDVI_bin.jpg", Mth);
}

void Run_NDWI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalNDWI(M[1],M[3]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("NDWI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_CIVIDIS);
    imwrite("NDWI_bin.jpg", Mth);
}

void Run_NWI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalNWI(M[0], M[3],M[4],M[6]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("NWI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_CIVIDIS);
    imwrite("NWI_bin.jpg", Mth);
}

void Run_DBI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalDBI(M[3], M[6]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("DBI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_CIVIDIS);
    imwrite("DBI_bin.jpg", Mth);
}

void Run_NDBI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalNDBI(M[3], M[6]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("NDBI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_CIVIDIS);
    imwrite("NDBI_bin.jpg", Mth);
}

void Run_IBI(Mat M[7]) {
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mnew = CalIBI(M[1],M[2],M[3], M[6]);
    Mat Mfeature = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    Mfeature = LineTransform(Mnew, 0, 255);
    Mfeature.convertTo(Mfeature, CV_8UC1);
    imshow("image", Mfeature);
    waitKey(0);
    imwrite("IBI.jpg", Mfeature);
    StateMethod2(Mfeature);
    Mat Mth = Mat::zeros(M[1].rows, M[1].cols, CV_8UC1);
    applyColorMap(Mfeature, Mth, COLORMAP_CIVIDIS);
    imwrite("IBI_bin.jpg", Mth);
}











int main()
{
    char imageChannel1[] = "./data/tm1.tif";
    char imageChannel2[] = "./data/tm2.tif";
    char imageChannel3[] = "./data/tm3.tif";
    char imageChannel4[] = "./data/tm4.tif";
    char imageChannel5[] = "./data/tm5.tif";
    char imageChannel6[] = "./data/tm6.tif";
    char imageChannel7[] = "./data/tm7.tif";

    Mat M[7];

    ReadImage(imageChannel1, M[0]);
    ReadImage(imageChannel2, M[1]);
    ReadImage(imageChannel3, M[2]);
    ReadImage(imageChannel4, M[3]);
    ReadImage(imageChannel5, M[4]);
    ReadImage(imageChannel6, M[5]);
    ReadImage(imageChannel7, M[6]);
    //ԭͼ����uchar16���ͣ�����ָ�����и����������Ҫת����ʽ
    for (int i = 0; i < 7; i++) {
        M[i].convertTo(M[i], CV_64FC1);
    }
    Run_RVI(M);
    Run_NDVI(M);
    Run_SAVI(M);
    Run_SAVI2(M);
    Run_TGDVI(M);
    Run_NDWI(M);
    Run_NWI(M);
    Run_DBI(M);
    Run_NDBI(M);
    Run_IBI(M);
}

