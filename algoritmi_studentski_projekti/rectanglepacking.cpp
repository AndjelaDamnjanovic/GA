#include "rectanglepacking.h"
#include "pomocnefunkcije.h"

#include <algorithm>
#include <QPainterPath>
#include<QSet>
#include <iostream>
#include <fstream>
#include <cmath>

#include "algoritambaza.h"

Rectangle::Rectangle(QPointF &topLeft, QPointF &bottomRight)
    : QRectF(topLeft, bottomRight)
{ }


/* Konstrukcija pravougaonika od parametara */
Rectangle::Rectangle(float x, float y, float w, float h)
    : QRectF(x, y, w, h)
{ }

rectanglePacking::rectanglePacking(QWidget *pCrtanje,
                                   int pauzaKoraka,
                                   const bool &naivni,
                                   std::string imeDatoteke,
                                   int brojPravougaonika,
                                   int poluprecnik)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni)
{
    if (imeDatoteke == ""){
        _rectangles = generateRandomRectangles(brojPravougaonika);
    }
    else
        _rectangles = readFromFile(imeDatoteke);
    _poluprecnik = poluprecnik;
    int xMin, xMax, yMin, yMax;
    initializeXsAndYs(&xMax, &xMin, &yMax, &yMin);
    _screenCenter = QPointF((xMax + xMin)/2, (yMax + yMin)/2);
};

void rectanglePacking::pokreniAlgoritam(){

    // povrsina najboljeg pakovanja
    int Pbest = 0;

    //AlgoritamBaza_updateCanvasAndBlock()

    while (_pCrtanje && Pbest == 0)
    {
        // ovde napisi sta treba da se radi pre crtanja
        //sortira OPADAJUCE
        //AlgoritamBaza_updateCanvasAndBlock()
        std::sort(_rectangles.begin(), _rectangles.end(), [](const Rectangle* a, const Rectangle* b)
                  { return a->height() * a->width()  > b->height() * b->width();});
        // procedura za pakovanje
        //QVector<Rectangle*> pbest = pack(_rectangles, _poluprecnik);
        // najbolje pakovanje
        QVector<Rectangle*> res;

        _circleCenter = QPointF(0, 0);

        // napravi sve moguce konfiguracuje E - onih 5 koje se pominju
        //za svaku konfiguraciju iz E:
        //napravi pakovanje i inicijalizuj granicu
        QVector <QVector <QPointF>> borders;
        QVector <QVector <Rectangle*> > configs = makeInitialConfigs(_rectangles, &borders);
        _configurations = configs;

        // dok ne dodjes do kraja liste nadji dozvoljeno mesto
        int n = configs.size();
        for(int i=0; i<n; i++){
            for(int j = configs[i].size(); j < _rectangles.size(); j++){
                // ako postoji, dodaj pravougaonik i updateuj granicu
                findFeasiblePlacement(_rectangles[j], &configs[i], &borders[i]);
                _configurations = configs;
                for (auto config : _configurations){
                    _res = shiftOnScreen(config);
                    AlgoritamBaza_updateCanvasAndBlock()
                }
            }
            //_res = shiftOnScreen(res);
            //AlgoritamBaza_updateCanvasAndBlock()
            // ako je p bolje od pbest onda je pbest = p, inace nista (po P ili br pravougaonika)
            qreal P = calcP(&configs[i]);
            if (P >= Pbest){
                Pbest = P;
                res = configs[i];
            }
        }
        _res = shiftOnScreen(res);
        AlgoritamBaza_updateCanvasAndBlock()
    }
    emit animacijaZavrsila();

}
void rectanglePacking::crtajAlgoritam(QPainter *painter) const{
    if (!painter) return;

    QPen p = painter->pen();
    p.setColor(Qt::red);
    p.setWidth(2);
    p.setCapStyle(Qt::RoundCap);

    painter->setPen(p);
    painter->drawEllipse(_screenCenter, _poluprecnik, _poluprecnik);

    p.setColor(Qt::blue);
    p.setWidth(5);
    painter->setPen(p);
    //nacrtaj pravougaonike
        for(auto rect : _res)
        {
            painter->drawRect(*rect);
            //painter->fillRect(rect);
        }
};


void rectanglePacking::pokreniNaivniAlgoritam(){
    emit animacijaZavrsila();
}
void rectanglePacking::crtajNaivniAlgoritam(QPainter *painter) const{
    if (!painter) return;
}

QVector<Rectangle*> rectanglePacking::readFromFile(const std::string imeDatoteke) const{
    QVector<Rectangle*> rects;
    std::ifstream inputFile(imeDatoteke);
    qreal x1, y1, x2, y2;
    while(inputFile >> x1 >> y1 >> x2 >> y2){
        QPointF tl = QPointF(x1, y1);
        QPointF br = QPointF(x2, y2);
        Rectangle* rect = new Rectangle(tl, br);
        rects.append(rect);
    }
    return rects;
}

QVector<Rectangle *> rectanglePacking::pack(QVector<Rectangle *> rectangles, int poluprecnik)
{
    // najbolje pakovanje
    QVector<Rectangle*> res;
    // povrsina najboljeg pakovanja
    int Pbest = 0;

    int xMin, xMax, yMin, yMax;
    initializeXsAndYs(&xMax, &xMin, &yMax, &yMin);

    _circleCenter = QPointF(0, 0);
    _screenCenter = QPointF((xMax + xMin)/2, (yMax + yMin)/2);
    _poluprecnik = poluprecnik;
    // napravi sve moguce konfiguracuje E - onih 5 koje se pominju
    //za svaku konfiguraciju iz E:
            //napravi pakovanje i inicijalizuj granicu
    QVector <QVector <QPointF>> borders;
    QVector <QVector <Rectangle*> > configs = makeInitialConfigs(_rectangles, &borders);

            // dok ne dodjes do kraja liste nadji dozvoljeno mesto
    int n = configs.size();
    for(int i=0; i<n; i++){
        for(int j = configs[i].size(); j < _rectangles.size(); j++){
            // ako postoji, dodaj pravougaonik i updateuj granicu
            findFeasiblePlacement(_rectangles[j], &configs[i], &borders[i]);
        }
        // ako je p bolje od pbest onda je pbest = p, inace nista (po P ili br pravougaonika)
        qreal P = calcP(&configs[i]);
        if (P >= Pbest){
            Pbest = P;
            res = configs[i];
        }
    }

    // vrati pbest
    return res;
}

void rectanglePacking::initializeXsAndYs(int * xMax, int *xMin, int *yMax, int *yMin) const
{
    if (_pCrtanje)
    {
        *xMax = _pCrtanje->width() - DRAWING_BORDER;
        *yMax = _pCrtanje->height() - DRAWING_BORDER;
    }
    else
    {
        *xMax = CANVAS_WIDTH;
        *yMax = CANVAS_HEIGHT;
    }

    *xMin = DRAWING_BORDER;
    *yMin = DRAWING_BORDER;
}

QVector<QVector<Rectangle *> > rectanglePacking::makeInitialConfigs(QVector<Rectangle *> rectangles, QVector <QVector <QPointF> > *borders) const
{
    QVector <QVector <Rectangle*>> configs;
    // konfiguracije idu redom: C, R, T, RS, RA
    // TREBA PROVERITI DA LI PRAVOUGAONIK UOPŠTE MOŽE DA STANE U KRUG NA POČETKU!!!!!!!



    while(intersectsCircle(rectangles.at(0)) && !rectangles.empty()){
        rectangles.removeFirst();
    }


    QVector<Rectangle*> config;

    if(rectangles.size() > 0){
        //config = makeConfigC(rectangles.at(0), borders);
        //configs.append(config);
        //config = makeConfigR(rectangles.at(0), borders);
        //configs.append(config);
        //config = makeConfigT(rectangles.at(0), borders);
        //configs.append(config);
    }
    if(rectangles.size() >= 2){
        config = makeConfigRA(rectangles.at(0), rectangles.at(1), borders);
        configs.append(config);
        //config = makeConfigRS(rectangles.at(0), rectangles.at(1), borders);
        //configs.append(config);
    }
    return configs;
}

QVector<Rectangle *> rectanglePacking::generateRandomRectangles(int numRectangles) const
{
    srand(static_cast<unsigned>(time(nullptr)));
    int xMin, xMax, yMin, yMax;
    initializeXsAndYs(&xMax, &xMin, &yMax, &yMin);
    QVector<Rectangle*> rects;

    int xDiff = xMax-xMin;
    int yDiff = yMax-yMin;
    Rectangle* rect;
    for(int i = 0; i < numRectangles; i++){
        QPointF fst = QPointF(xMin + rand()%xDiff, yMin + rand()%yDiff);
        QPointF snd =  QPointF(xMin + rand()%xDiff, yMin + rand()%yDiff);

        if (fst.x() < snd.x()){
            if (fst.y() < snd.y()){
                rect = new Rectangle(fst, snd);
            }else{
                QPointF left = QPointF(fst.x(), snd.y());
                QPointF right = QPointF(snd.x(), fst.y());
                rect = new Rectangle(left, right);
            }
        }else{
            if(snd.y() > fst.y())
                rect = new Rectangle(snd, fst);
            else{
                QPointF left = QPointF(snd.x(), fst.y());
                QPointF right = QPointF(fst.x(), snd.y());
                rect = new Rectangle(left, right);
            }
        }
        rects.append(rect);
    }

    return rects;
}
QVector <Rectangle*> rectanglePacking::makeConfigC(Rectangle* rect, QVector <QVector <QPointF> > *borders) const{
    QVector <Rectangle*> config;
    QPointF center = _circleCenter;

    // premesti centar pravougaonika u centar kruga
    QPointF rectCenter = rect->center();
    float xDiff = rectCenter.x() - center.x();
    float yDiff = rectCenter.y() - center.y();

    qreal x1, y1, x2, y2;
    rect->getCoords(&x1, &y1, &x2, &y2);
    QPointF tl = QPointF(x1 - xDiff, y1 - yDiff);
    QPointF br = QPointF(x2 - xDiff, y2 - yDiff);
    Rectangle * resRect = new Rectangle(tl, br);
    config.append(resRect);

    resRect->getCoords(&x1, &y1, &x2, &y2);
    QVector <QPointF> border;
    border.append(QPointF(x2, y1));
    border.append(QPointF(x1, y1));
    border.append(QPointF(x1, y2));
    border.append(QPointF(x2, y2));
    borders->append(border);
    return config;
};
QVector <Rectangle*> rectanglePacking::makeConfigR(Rectangle* rect, QVector <QVector <QPointF> >* borders) const{
    QVector <Rectangle*> config;
    QPointF center = _circleCenter;
    qreal L = rect->width();
    qreal W = rect->height();

    //siftuj pravougaonik na desno
    if((std::pow(_poluprecnik, 2) - std::pow(W, 2)/4) >= 0){
        QPointF rectCenter = rect->center();
        QPointF rectCenterNew = QPointF(std::sqrt(std::pow(_poluprecnik, 2) - std::pow(W, 2)/4) - L/2, 0);

        float xDiff = rectCenterNew.x() - rectCenter.x();
        float yDiff = rectCenterNew.y() - rectCenter.y();

        qreal x1, y1, x2, y2;
        rect->getCoords(&x1, &y1, &x2, &y2);
        QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
        QPointF br = QPointF(x2 + xDiff, y2 + yDiff);

        Rectangle * resRect = new Rectangle(tl, br);
        resRect->getCoords(&x1, &y1, &x2, &y2);
        QVector <QPointF> border;
        border.append(QPointF(x2, y1));
        border.append(QPointF(x1, y1));
        border.append(QPointF(x1, y2));
        border.append(QPointF(x2, y2));
        borders->append(border);

        config.append(resRect);
    }
    return config;
};
QVector <Rectangle*> rectanglePacking::makeConfigT(Rectangle* rect, QVector <QVector <QPointF> >* borders) const{
    QVector <Rectangle*> config;
    QPointF center = _circleCenter;

    qreal L = rect->width();
    qreal W = rect->height();

    // siftuj pravougaonik na gore
    if(std::pow(_poluprecnik, 2) - std::pow(L, 2)/4 >= 0){
        QPointF rectCenter = rect->center();
        QPointF rectCenterNew = QPointF(0, std::sqrt(std::pow(_poluprecnik, 2) - std::pow(L, 2)/4) - W/2);

        float xDiff = rectCenterNew.x() - rectCenter.x();
        float yDiff = rectCenterNew.y() - rectCenter.y();

        qreal x1, y1, x2, y2;
        rect->getCoords(&x1, &y1, &x2, &y2);

        QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
        QPointF br = QPointF(x2 + xDiff, y2 + yDiff);

        Rectangle * resRect = new Rectangle(tl, br);
        resRect->getCoords(&x1, &y1, &x2, &y2);

        QVector <QPointF> border;
        border.append(QPointF(x2, y1));
        border.append(QPointF(x1, y1));
        border.append(QPointF(x1, y2));
        border.append(QPointF(x2, y2));
        borders->append(border);

        config.append(resRect);
    }

    return config;
};
QVector <Rectangle*> rectanglePacking::makeConfigRS(Rectangle* fst, Rectangle* snd, QVector <QVector <QPointF> >* borders) const{
    QVector <QVector <QPointF>>* borders2 = new QVector <QVector <QPointF>>;
    QVector <Rectangle*> config = makeConfigR(fst, borders2);
    //config.append(fst);
    QPointF center = _circleCenter;
    Rectangle* rect = config.at(0);
    if(intersectsCircle(rect)){
        //return config;
    }
    qreal x1, y1, x2, y2;
    qreal x11, y11, x12, y12;
    rect->getCoords(&x11, &y11, &x12, &y12);

    qreal L1 = fst->width();
    qreal W1 = fst->height();

    qreal L2 = snd->width();
    qreal W2 = snd->height();

    QPointF fstCenter = rect->center();

    if(std::pow(_poluprecnik, 2) - std::pow(fstCenter.x() - L1/2 - L2,2) >= 0){
        QPointF rectCenter = snd->center();

        //QPoint rectCenter = QPoint(0,0);
        QPointF rectCenterNew = QPointF(fstCenter.x() - L1/2 - L2/2,
                                        std::sqrt(std::pow(_poluprecnik, 2) - std::pow(fstCenter.x() - L1/2 - L2, 2)) - W2/2);


        qreal xDiff = rectCenterNew.x() - rectCenter.x();
        qreal yDiff = rectCenterNew.y() - rectCenter.y();

        snd->getCoords(&x1, &y1, &x2, &y2);
        QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
        QPointF br = QPointF(x2 + xDiff, y2 + yDiff);

        Rectangle * resRect = new Rectangle(tl, br);
        if(intersectsCircle(rect)){
            //return config;
        }
        // proveri da li staje u krug
        //config.append(snd);
        if(rectCenterNew.y() <= (W1 + W2)/2){
            config.append(resRect);
            QVector <QPointF> border;
            resRect->getCoords(&x1, &y1, &x2, &y2);
            //PROVERI DA SE TEMENA NE POKLAPAJU!!!!
            border.append(QPointF(x12, y11));
            border.append(QPointF(x11, y11));
            border.append(QPointF(x2, y1));
            border.append(QPointF(x1, y1));
            border.append(QPointF(x1, y2));
            border.append(QPointF(x2, y2));
            border.append(QPointF(x11, y12));
            border.append(QPointF(x12, y12));
            borders->append(border);
        }else{
            //return config;
        }
    }else{
       // return config;
    }
    return config;
};
QVector <Rectangle*> rectanglePacking::makeConfigRA(Rectangle* fst, Rectangle* snd, QVector <QVector <QPointF> >* borders) const{
    QVector <Rectangle*> config;
    //config.append(&fst);
    QPointF center = _circleCenter;

    qreal L1 = fst->width();
    qreal W1 = fst->height();
    qreal L2 = snd->width();
    qreal W2 = snd->height();

    QPointF rectCenterSnd = snd->center();
    QPointF rectCenterFst = fst->center();

    QPointF fstCenter = QPointF((W1 + W2)*std::sqrt((4*std::pow(_poluprecnik, 2)/(std::pow(L1 - L2, 2) + std::pow(W1 + W2, 2)) - 1))/2 - L2/2, W2/2 - (L1 - L2)*std::sqrt((4*std::pow(_poluprecnik, 2)/(std::pow(L1 -L2, 2) + std::pow(W1 + W2, 2)))-1)/2);
    QPointF sndCenter = QPointF(fstCenter.x() - L1/2 + L2/2, fstCenter.y() - W1/2 -W2/2);

    if((std::pow(L1/2 -L2/2, 2) + std::pow(W1/2 + W2/2, 2)) <= std::pow(_poluprecnik, 2)){
        float xDiff = fstCenter.x() - rectCenterFst.x();
        float yDiff = fstCenter.y() - rectCenterFst.y();

        qreal x1, y1, x2, y2;
        qreal x11, y11, x12, y12;
        fst->getCoords(&x11, &y11, &x12, &y12);

        QPointF tl = QPointF(x11 + xDiff, y11 + yDiff);
        QPointF br = QPointF(x12 + xDiff, y12 + yDiff);

        Rectangle * resRect = new Rectangle(tl, br);
        std::cout<<"RA konfig prvi pravougaonik: "<<std::to_string(fstCenter.x())<<" "<<std::to_string(fstCenter.y())<<std::endl;
        if(intersectsCircle(resRect)){
            //return config;
        }
        resRect->getCoords(&x11, &y11, &x12, &y12);

        config.append(resRect);

        xDiff = sndCenter.x() - rectCenterSnd.x();
        yDiff = sndCenter.y() - rectCenterSnd.y();

        snd->getCoords(&x1, &y1, &x2, &y2);
        tl = QPointF(x1 + xDiff, y1 + yDiff);
        br = QPointF(x2 + xDiff, y2 + yDiff);

        resRect = new Rectangle(tl, br);
        if(intersectsCircle(resRect)){
            //return config;
        }

        resRect->getCoords(&x1, &y1, &x2, &y2);

        config.append(resRect);

        QVector <QPointF> border;

        border.append(QPointF(x12, y11));
        border.append(QPointF(x2, y2));
        border.append(QPointF(x2, y1));
        border.append(QPointF(x1, y1));
        border.append(QPointF(x11, y12));
        border.append(QPointF(x12, y12));
        borders->append(border);
    }/*else{
        return config;
    }*/
    return config;
}

void rectanglePacking::findFeasiblePlacement(Rectangle *rect, QVector<Rectangle *> *config, QVector<QPointF> *border) const
    {
    QVector<QVector <QPointF>> possibleBorders;
    QVector<Rectangle*> possiblePositions;
    qreal minDist = INT_MAX;

    for(int i=0; i<border->size() - 2; i++){
        QVector<QPointF> newBorder = *border;
        QPointF A = border->at(i);
        QPointF B = border->at(i+1);
        QPointF C = border->at(i+2);

        QPointF BC = QPointF(C.x() - B.x(), C.y() - B.y());
        QPointF CB = QPointF(B.x() - C.x(), B.y() - C.y());
        QPointF BA = QPointF(A.x() - B.x(), A.y() - B.y());
        QPointF AB = QPointF(B.x() - A.x(), B.y() - A.y());
        QPointF centerNew;


        qreal BCnorm = std::sqrt(std::pow(BC.x(), 2) + std::pow(BC.y(), 2));
        qreal CBnorm = std::sqrt(std::pow(CB.x(), 2) + std::pow(CB.y(), 2));
        qreal ABnorm = std::sqrt(std::pow(AB.x(), 2) + std::pow(AB.y(), 2));
        qreal BAnorm = std::sqrt(std::pow(BA.x(), 2) + std::pow(BA.y(), 2));

        qreal L = rect->width();
        qreal W = rect->height();

        qreal crossProduct = (B.x() - A.x())*(C.y() - B.y()) - (B.y() - A.y())*(C.x() - B.x());
        if ( crossProduct > 0 ){
            //C je desno od AB
            //IMA GRESKE U GRANICI BAR U SLUCAJU 1
            processCaseA(&possibleBorders, &possiblePositions, &newBorder, rect, L, W, &BA, &BC, BAnorm, BCnorm, &A, &B, &C, i);

        }else{
                // slucaj b)
                processCaseB(&possibleBorders, &possiblePositions, &newBorder, rect, L, W, &BA, &CB, BAnorm, CBnorm, &A, &B, &C, i);
                // slucaj c)
                processCaseC(&possibleBorders, &possiblePositions, &newBorder, rect, L, W, &BC, &AB, BCnorm, ABnorm, &A, &B, &C, i);
                //processCaseD(&possibleBorders, &possiblePositions, &newBorder, rect, L, W, &BA, &CB, BAnorm, CBnorm, &A, &B, &C, i);
        }
    }

    // ako ima dozvoljenih pozicija, gledamo koja je najbliza centru i nju dodajemo
    int n = possiblePositions.size();

    if(n > 0){
        int index = findNearestToCenter(&possiblePositions);
        *border = possibleBorders[index];
        config->append(possiblePositions[index]);
    }
    //return;
    };

bool rectanglePacking::intersectsCircle(Rectangle * resRect) const
{
    qreal x1, x2, y1, y2;
    resRect->getCoords(&x1, &y1, &x2, &y2);
    //proveri sece li se neko teme s krugom -- uzmi i iyracunaj rastojanje svakog temena od (0,0) i ako je vise od 287 sece se
    qreal tlNorm = std::sqrt(std::pow(x1, 2) + std::pow(y1, 2));
    qreal brNorm = std::sqrt(std::pow(x2, 2) + std::pow(y2, 2));
    qreal trNorm = std::sqrt(std::pow(x2, 2) + std::pow(y1, 2));
    qreal blNorm = std::sqrt(std::pow(x1, 2) + std::pow(y2, 2));
    if((tlNorm > 287) || (blNorm > 287) || (trNorm > 287) || (brNorm > 287)){
        return true;
    }
    return false;
}

void rectanglePacking::processCaseA(QVector<QVector<QPointF> > *possibleBorders, QVector<Rectangle *> *possiblePositions, QVector<QPointF> * newBorder, Rectangle* rect, qreal L, qreal W, QPointF *BA, QPointF *BC, qreal BAnorm, qreal BCnorm, QPointF* A, QPointF* B, QPointF* C, int i) const
{
    std::cout<<"U process A sam!"<<std::endl;
    QPointF centerNew;
    //Ako je AB horizontalna
    if ( A->y() == B->y() )
        centerNew = QPointF(B->x() + (L/(2*BAnorm))* BA->x() + (W/(2*BCnorm))*BC->x(), B->y() + (L/(2*BAnorm))* BA->y() + (W/(2*BCnorm))*BC->y());
    else
        centerNew = QPointF(B->x() + (W/(2*BAnorm))* BA->x() + (L/(2*BCnorm))*BC->x(), B->y() + (W/(2*BAnorm))* BA->y() + (L/(2*BCnorm))*BC->y());
    QPointF center = rect->center();

    qreal xDiff = centerNew.x() - center.x();
    qreal yDiff = centerNew.y() - center.y();

    qreal x1, x2, y1, y2;
    rect->getCoords(&x1, &y1, &x2, &y2);

    QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
    QPointF br = QPointF(x2 + xDiff, y2 + yDiff);

    Rectangle* resRect = new Rectangle(tl, br);
    std::cout<<"Res rect u A:"<< std::to_string(resRect->center().x())<<" "<<std::to_string(resRect->center().y());

    if (intersectsCircle(resRect)){
        return;
    }else{
    newBorder->removeAt(i+1);
    // u zavisnosti od polozaja tacaka A i C, moze se videti tacno koje tacke treba dodati u granicu
    if (C->x() < A->x()){
        if (C->y() < A->y()){
                newBorder->insert(i+1, br);
                newBorder->insert(i+2, QPointF(br.x(), tl.y()));
                newBorder->insert(i+3, tl);
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }else{
                newBorder->insert(i+1, QPointF(br.x(), tl.y()));
                newBorder->insert(i+2, tl);
                newBorder->insert(i+3, QPointF(tl.x(), br.y()));
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }
    }else{
        if(A->y() < C->y()){
                newBorder->insert(i+1, tl);
                newBorder->insert(i+2, QPointF(tl.x(), br.y()));
                newBorder->insert(i+3, br);
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }else{
                newBorder->insert(i+1, QPointF(tl.x(), br.y()));
                newBorder->insert(i+2, br);
                newBorder->insert(i+3, QPointF(br.x(), tl.y()));
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }
    }
    }
};

void rectanglePacking::processCaseB(QVector<QVector<QPointF> > *possibleBorders, QVector<Rectangle *> *possiblePositions, QVector<QPointF> * newBorder, Rectangle* rect, qreal L, qreal W, QPointF *BA, QPointF *CB, qreal BAnorm, qreal CBnorm, QPointF* A, QPointF* B, QPointF* C, int i) const
{
    QPointF centerNew;
    std::cout<<"U process B sam!"<<std::endl;
    //Ako je AB horizontalna
    if ( A->y() == B->y() )
        centerNew = QPointF(B->x() + (L/(2*BAnorm))* BA->x() + (W/(2*CBnorm))*CB->x(), B->y() + (L/(2*BAnorm))* BA->y() + (W/(2*CBnorm))*CB->y());
    else
        centerNew = QPointF(B->x() + (W/(2*BAnorm))* BA->x() + (L/(2*CBnorm))*CB->x(), B->y() + (W/(2*BAnorm))* BA->y() + (L/(2*CBnorm))*CB->y());
    QPointF center = rect->center();

    qreal xDiff = centerNew.x() - center.x();
    qreal yDiff = centerNew.y() - center.y();

    qreal x1, x2, y1, y2;
    rect->getCoords(&x1, &y1, &x2, &y2);

    QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
    QPointF br = QPointF(x2 + xDiff, y2 + yDiff);
    Rectangle* resRect = new Rectangle(tl, br);

    if (intersectsCircle(resRect)){
        processCaseD(possibleBorders, possiblePositions, newBorder, resRect, L, W, BA, CB, BAnorm, CBnorm, A, B, C, i);
        return;
    }else{
        newBorder->removeAt(i+1);
        // u zavisnosti od polozaja tacaka A i C, moze se videti tacno koje tacke treba dodati u granicu

        if (C->x() < A->x()){
            if (C->y() < A->y()){
                    newBorder->insert(i+1, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+2, br);
                    newBorder->insert(i+3, QPointF(br.x(), tl.y()));
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
            }else{
                    newBorder->insert(i+1, br);
                    newBorder->insert(i+2, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+3, tl);
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
            }
        }else{
            if(A->y() < C->y()){
                    newBorder->insert(i+1, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+2, tl);
                    newBorder->insert(i+3, QPointF(tl.x(), br.y()));
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
            }else{
                    newBorder->insert(i+1, tl);
                    newBorder->insert(i+2, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+3, br);
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
            }
        }
    }
};

void rectanglePacking::processCaseC(QVector<QVector<QPointF> > *possibleBorders, QVector<Rectangle *> *possiblePositions, QVector<QPointF> * newBorder, Rectangle* rect, qreal L, qreal W, QPointF *BC, QPointF *AB, qreal BCnorm, qreal ABnorm, QPointF* A, QPointF* B, QPointF* C, int i) const
{
    QPointF centerNew;
    std::cout<<"U process C sam!"<<std::endl;
    //Ako je AB horizontalna
    if ( A->y() == B->y() )
        centerNew = QPointF(B->x() + (W/(2*BCnorm))* BC->x() + (L/(2*ABnorm))*AB->x(), B->y() + (W/(2*BCnorm))* BC->y() + (L/(2*ABnorm))*AB->y());
    else
        centerNew = QPointF(B->x() + (L/(2*BCnorm))* BC->x() + (W/(2*ABnorm))*AB->x(), B->y() + (L/(2*BCnorm))* BC->y() + (W/(2*ABnorm))*AB->y());
    QPointF center = rect->center();

    qreal xDiff = centerNew.x() - center.x();
    qreal yDiff = centerNew.y() - center.y();

    qreal x1, x2, y1, y2;
    rect->getCoords(&x1, &y1, &x2, &y2);

    QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
    QPointF br = QPointF(x2 + xDiff, y2 + yDiff);
    Rectangle* resRect = new Rectangle(tl, br);
    std::cout<<"Res rect u C:"<< std::to_string(resRect->topLeft().x())<<" "<<std::to_string(resRect->topLeft().y());

    if (intersectsCircle(resRect)){
        processCaseE(possibleBorders, possiblePositions, newBorder, resRect, L, W, BC, AB, BCnorm, ABnorm, A, B, C, i);
        return;
    }
    newBorder->removeAt(i+1);
    // u zavisnosti od polozaja tacaka A i C, moze se videti tacno koje tacke treba dodati u granicu
    if (C->x() < A->x()){
        if (C->y() < A->y()){
                newBorder->insert(i+1, QPointF(br.x(), tl.y()));
                newBorder->insert(i+2, tl);
                newBorder->insert(i+3, QPointF(tl.x(), br.y()));
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }else{
                newBorder->insert(i+1, tl);
                newBorder->insert(i+2, QPointF(tl.x(), br.y()));
                newBorder->insert(i+3, br);
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }
    }else{
        if(A->y() < C->y()){
                newBorder->insert(i+1, QPointF(tl.x(), br.y()));
                newBorder->insert(i+2, br);
                newBorder->insert(i+3, QPointF(br.x(), tl.y()));
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }else{
                newBorder->insert(i+1, br);
                newBorder->insert(i+2, QPointF(br.x(), tl.y()));
                newBorder->insert(i+3, tl);
                possibleBorders->append(*newBorder);
                possiblePositions->append(resRect);
        }
    }
};

void rectanglePacking::processCaseD(QVector<QVector<QPointF> > *possibleBorders, QVector<Rectangle *> *possiblePositions, QVector<QPointF> * newBorder, Rectangle* rect, qreal L, qreal W, QPointF *BA, QPointF *CB, qreal BAnorm, qreal CBnorm, QPointF* A, QPointF* B, QPointF* C, int i) const
{
    qreal lambda;
    std::cout<<"U process D sam!"<<std::endl;
    QPointF centerNew;
    //Ako je AB horizontalna
    if ( A->y() == B->y() ){
        lambda = calculateLambdaHorizontal(BA, BAnorm, CB, CBnorm, W, B);
        if (lambda != INT_MAX)
                centerNew = QPointF(rect->center() + lambda*(*BA)/BAnorm);
    }else{
        lambda = calculateLambdaVertical(BA, BAnorm, CB, CBnorm, L, B);
        if (lambda != INT_MAX)
                centerNew = QPointF(rect->center() + lambda*(*BA)/BAnorm);
    }


    if(centerNew != QPointF(0,0)){
        QPointF center = rect->center();
        qreal xDiff = centerNew.x() - center.x();
        qreal yDiff = centerNew.y() - center.y();

        qreal x1, x2, y1, y2;
        rect->getCoords(&x1, &y1, &x2, &y2);

        QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
        QPointF br = QPointF(x2 + xDiff, y2 + yDiff);
        Rectangle* resRect = new Rectangle(tl, br);

        std::cout<<"Res rect u D:"<< std::to_string(resRect->center().x())<<" "<<std::to_string(resRect->center().y());
        if (intersectsCircle(resRect)){
                return;
        }

        // u zavisnosti od polozaja tacaka A i C, moze se videti tacno koje tacke treba dodati u granicu
        if (C->x() < A->x()){
                if (C->y() < A->y()){
                    newBorder->insert(i+1, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+2, br);
                    newBorder->insert(i+3, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+4, tl);
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }else{
                    newBorder->insert(i+1, br);
                    newBorder->insert(i+2, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+3, tl);
                    newBorder->insert(i+4, QPointF(tl.x(), br.y()));
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }
        }else{
                if(A->y() < C->y()){
                    newBorder->insert(i+1, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+2, tl);
                    newBorder->insert(i+3, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+4, br);
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }else{
                    newBorder->insert(i+1, tl);
                    newBorder->insert(i+2, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+3, br);
                    newBorder->insert(i+4, QPointF(br.x(), tl.y()));
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }
        }
    }

};

qreal rectanglePacking::calculateLambdaHorizontal(QPointF * BA, qreal BAnorm, QPointF *CB, qreal CBnorm, qreal W, QPointF* B) const
{
    qreal lambda;
    QPointF coefficient = (*B) + W*(*CB)/CBnorm;
    QPointF lambdaMultiplier = (*BA)/BAnorm;

    qreal x = coefficient.x();
    qreal y = coefficient.y();

    //lambdaMultiplier/=denominator;

    qreal a = std::pow(lambdaMultiplier.x(), 2) + std::pow(lambdaMultiplier.y(), 2);
    qreal b = lambdaMultiplier.x() * x + lambdaMultiplier.y() * y;
    qreal c = std::pow(x, 2) + std::pow(y, 2) - std::pow(_poluprecnik, 2);

    //resi kvadratnu jednacinu a*(lambda^2) + 2*b*lambda + c = 0
    if(std::pow(b, 2) - a*c > 0){
        //ima resenja
        qreal lambda1 = (std::sqrt(std::pow(b, 2) - a*c) -b)/a;
        qreal lambda2 = (-std::sqrt(std::pow(b,2) - a*c) - b)/a;

        if(std::abs(lambda1) < std::abs(lambda2)){
                lambda = lambda1;
        }else{
                lambda = lambda2;
        }

        //proveri da li lambda zadovoljava uslove i ako da, vrati ga
        if(lambda < 0){
                if(std::abs(lambda) <= W){
                   return lambda;
                }else{
                   return INT_MAX;
                }
        }else{
                if(lambda <= CBnorm){
                   return lambda;
                }else{
                   return INT_MAX;
                }
        }
    }
    return INT_MAX;
};

qreal rectanglePacking::calculateLambdaVertical(QPointF * BA, qreal BAnorm, QPointF * CB, qreal CBnorm, qreal L, QPointF* B) const
{
    qreal lambda;
    QPointF coefficient = (*B) + L*(*CB)/CBnorm;
    QPointF lambdaMultiplier = (*BA)/BAnorm;

    qreal x = coefficient.x();
    qreal y = coefficient.y();

    //lambdaMultiplier/=denominator;

    qreal a = std::pow(lambdaMultiplier.x(), 2) + std::pow(lambdaMultiplier.y(), 2);
    qreal b = lambdaMultiplier.x() * x + lambdaMultiplier.y() * y;
    qreal c = std::pow(x, 2) + std::pow(y, 2) - std::pow(_poluprecnik, 2);

    //resi kvadratnu jednacinu a*(lambda^2) + 2*b*lambda + c = 0
    if(std::pow(b, 2) - a*c > 0){
        //ima resenja
        qreal lambda1 = (std::sqrt(std::pow(b, 2) - a*c) -b)/a;
        qreal lambda2 = (-std::sqrt(std::pow(b,2) - a*c) - b)/a;

        if(std::abs(lambda1) < std::abs(lambda2)){
                lambda = lambda1;
        }else{
                lambda = lambda2;
        }

        //proveri da li lambda zadovoljava uslove i ako da, vrati ga
        if(lambda < 0){
                if(std::abs(lambda) <= L){
                    return lambda;
                }else{
                    return INT_MAX;
                }
        }else{
                if(lambda <= CBnorm){
                    return lambda;
                }else{
                    return INT_MAX;
                }
            }
    }
    return INT_MAX;
};

int rectanglePacking::findNearestToCenter(QVector<Rectangle *> *candidates) const
{
    int n = candidates->size();
    int closestIndex = 0;
    qreal distMin = calcDistance((candidates->at(0))->center());
    for (int i = 1; i < n; i++){
            qreal currDist = calcDistance((candidates->at(i))->center());
        if (currDist < distMin)
            closestIndex = i;
    }
    return closestIndex;
};

qreal rectanglePacking::calcDistance(QPointF rectCenter) const
{
    return std::sqrt(std::pow(rectCenter.x(), 2) + std::pow(rectCenter.y(), 2));
};

qreal rectanglePacking::calcP(QVector<Rectangle *> *config) const
{
    qreal res = 0;
    for(auto rect : *config){
        res += rect->height() * rect->width();
    }
    return res;
};

QVector<Rectangle *> rectanglePacking::shiftOnScreen(QVector<Rectangle *> rectangles) const
{
    Rectangle* resRect;
    QVector<Rectangle*> res;
    for(auto rect : rectangles){
        qreal x1, x2, y1, y2;
        rect->getCoords(&x1, &y1, &x2, &y2);

        QPointF tl = QPointF(x1 + _screenCenter.x(), y1 + _screenCenter.y());
        QPointF br = QPointF(x2 + _screenCenter.x(), y2 + _screenCenter.y());
        resRect = new Rectangle(tl, br);
        res.append(resRect);
    }
    return res;
};

void rectanglePacking::processCaseE(QVector<QVector<QPointF> > *possibleBorders, QVector<Rectangle *> *possiblePositions, QVector<QPointF> * newBorder, Rectangle* rect, qreal L, qreal W, QPointF *BC, QPointF *AB, qreal BCnorm, qreal ABnorm, QPointF* A, QPointF* B, QPointF* C, int i) const
{
    std::cout<<"U process E sam!"<<std::endl;
    qreal lambda;

    QPointF centerNew;
    //Ako je AB horizontalna
    if ( A->y() == B->y() ){
            lambda = calculateLambdaHorizontal(BC, BCnorm, AB, ABnorm, L, B);
            if (lambda != INT_MAX)
                centerNew = QPointF(rect->center() + lambda*(*BC)/BCnorm);
    }else{
            lambda = calculateLambdaVertical(BC, BCnorm, AB, ABnorm, W, B);
            if (lambda != INT_MAX)
                centerNew = QPointF(rect->center() + lambda*(*BC)/BCnorm);
    }


    if(centerNew != QPointF(0,0)){
            QPointF center = rect->center();
            qreal xDiff = centerNew.x() - center.x();
            qreal yDiff = centerNew.y() - center.y();

            qreal x1, x2, y1, y2;
            rect->getCoords(&x1, &y1, &x2, &y2);

            QPointF tl = QPointF(x1 + xDiff, y1 + yDiff);
            QPointF br = QPointF(x2 + xDiff, y2 + yDiff);
            Rectangle* resRect = new Rectangle(tl, br);

            std::cout<<"Res rect u E:"<< std::to_string(resRect->topLeft().x())<<" "<<std::to_string(resRect->topLeft().y());
            if (intersectsCircle(resRect)){
                return;
            }

            // u zavisnosti od polozaja tacaka A i C, moze se videti tacno koje tacke treba dodati u granicu
            if (C->x() < A->x()){
                if (C->y() < A->y()){
                    newBorder->insert(i+2, br);
                    newBorder->insert(i+3, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+4, tl);
                    newBorder->insert(i+5, QPointF(tl.x(), br.y()));
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }else{
                    newBorder->insert(i+2, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+3, tl);
                    newBorder->insert(i+4, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+5, br);
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }
            }else{
                if(A->y() < C->y()){
                    newBorder->insert(i+2, tl);
                    newBorder->insert(i+3, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+4, br);
                    newBorder->insert(i+5, QPointF(br.x(), tl.y()));
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }else{
                    newBorder->insert(i+2, QPointF(tl.x(), br.y()));
                    newBorder->insert(i+3, br);
                    newBorder->insert(i+4, QPointF(br.x(), tl.y()));
                    newBorder->insert(i+5, tl);
                    possibleBorders->append(*newBorder);
                    possiblePositions->append(resRect);
                }
            }
    }

};
