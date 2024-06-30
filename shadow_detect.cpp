#include "head.h"
Mat HSVShadowTest(Mat Img) {
    Mat HSVMat;
    int nHeight = Img.rows;//��ȡͼ��ĸ�
    int nWidth = Img.cols;//��ȡͼ��Ŀ�
    int nChannels = Img.channels();//��ͨ��
    double* pImData = (double*)Img.data;
    Mat Mnew = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNew = (double*)Mnew.data;
    //������Ԫ
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double B = *(pImData + (nWidth * i + j) * nChannels);
            double G = *(pImData + (nWidth * i + j) * nChannels + 1);
            double R = *(pImData + (nWidth * i + j) * nChannels + 2);
            double V = (B + G + R) / 3;
            double theta = acos(((R - G) + (R - B)) / 2 / sqrt((R - G) * (R - G) + (R - B) * (G - B))) / 3.1415926 * 180;
            double H = theta;
            if (B > G) {
                H = 360 - theta;
            }
            double S = 1 - 3 * min(R, min(G, B)) / (R + G + B);
            double M = (S - V) / (H + S + V);
            *(pImDataNew + (nWidth * i + j)) = M;
        }
    }
    return Mnew;
}


Mat C1C2C3ShadowTest(Mat Img) {
    Mat HSVMat;
    int nHeight = Img.rows;//��ȡͼ��ĸ�
    int nWidth = Img.cols;//��ȡͼ��Ŀ�
    int nChannels = Img.channels();//��ȡͼ��ͨ����Ŀ BGR/��ͨ��
    double* pImData = (double*)Img.data;
    Mat Mnew = Mat::zeros(nHeight, nWidth, CV_8UC1);
    unsigned char* pImDataNew = Mnew.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double B = *(pImData + (nWidth * i + j) * nChannels);
            double G = *(pImData + (nWidth * i + j) * nChannels + 1);
            double R = *(pImData + (nWidth * i + j) * nChannels + 2);
            double C1 = atan(R / max(G, B));
            double C2 = atan(G / max(R, B));
            double C3 = atan(B / max(G, R));
            if (C3 > 0.4 && B < 30) {
                *(pImDataNew + (nWidth * i + j)) = 255;
            }
            else {
                *(pImDataNew + (nWidth * i + j)) = 0;
            }
        }
    }
    return Mnew;
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


//ʹ���ȸ�ʴ�����͵ĺ�����
Mat After_process(Mat Mnew) {
    Mat element1 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(3, 3));
    dilate(Mnew, Mnew, element1);
    Mat element2 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(4, 4));
    erode(Mnew, Mnew, element2);
    return Mnew;
}
Mat After_process2(Mat Mnew) {
    Mat element1 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(6, 6));
    dilate(Mnew, Mnew, element1);
    Mat element2 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(7, 7));
    erode(Mnew, Mnew, element2);
    return Mnew;
}

bool ReadImage(char imageName[], Mat& M) {
    M = imread(imageName, IMREAD_ANYCOLOR);
    if (M.empty())     // �ж��ļ��Ƿ�������     
    {
        fprintf(stderr, "Can not load image %s\n", imageName);
        waitKey(6000);  // �ȴ�6000 ms�󴰿��Զ��ر�   
        return false;
    }
}

void Run_HSV_detect() {
    char ShadowImageName[] = "./shadow_data/Color.bmp";
    Mat MShadow;
    ReadImage(ShadowImageName, MShadow);
    MShadow.convertTo(MShadow, CV_64F);
    Mat Mnew = Mat::zeros(MShadow.rows, MShadow.cols, CV_64FC1);
    Mnew = HSVShadowTest(MShadow);
    Mnew = LineTransform(Mnew, 0, 255);
    Mnew.convertTo(Mnew, CV_8UC1);
    Mat M2 = Mat::zeros(MShadow.rows, MShadow.cols, CV_64FC1);
    threshold(Mnew, M2, 127, 255, THRESH_OTSU);
    applyColorMap(M2, M2, COLORMAP_CIVIDIS);
    M2 = After_process(M2);
    imshow("image", M2);
    waitKey();
    imwrite("HSV_detect.jpg", M2);
}

void Run_C_detect() {
    char ShadowImageName[] = "./shadow_data/Color.bmp";
    Mat MShadow;
    ReadImage(ShadowImageName, MShadow);
    MShadow.convertTo(MShadow, CV_64F);
    Mat Mnew = Mat::zeros(MShadow.rows, MShadow.cols, CV_64FC1);
    Mnew = C1C2C3ShadowTest(MShadow);
    //Mnew = LineTransform(Mnew, 0, 255);
    Mnew.convertTo(Mnew, CV_8UC1);
    Mat M2 = Mat::zeros(MShadow.rows, MShadow.cols, CV_64FC1);
    threshold(Mnew, M2, 127, 255, THRESH_OTSU);
    applyColorMap(M2, M2, COLORMAP_CIVIDIS);
    M2 = After_process2(M2);
    imshow("image", M2);
    waitKey();
    imwrite("C_detect.jpg", M2);
}



int main() {
    Run_HSV_detect();
    Run_C_detect();
    return 0;
}

