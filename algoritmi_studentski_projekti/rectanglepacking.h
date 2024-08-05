#ifndef RECTANGLEPACKING_H
#define RECTANGLEPACKING_H

#include "algoritambaza.h"
#include <QVector>
#include <QSet>


struct Rectangle : QRectF {
    /* Konstruktori i destruktor strukture */
    Rectangle(QPointF &, QPointF &);
    Rectangle(float, float, float, float);
    virtual ~Rectangle() = default;
};

class rectanglePacking : public AlgoritamBaza
{
public:
    rectanglePacking(QWidget *pCrtanje,
                     int pauzaKoraka,
                     const bool &naivni = false,
                     std::string imeDatoteke = "",
                     int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA,
                     int poluprecnik = 50);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

    QVector<Rectangle*> generateRandomRectangles(int numRectangles = BROJ_SLUCAJNIH_OBJEKATA) const;
    QVector<Rectangle*> readFromFile(const std::string) const;
    QVector<Rectangle*> pack(QVector<Rectangle*>, int);
    void initializeXsAndYs(int*, int*, int*, int*) const;
    QVector <QVector <Rectangle*> > makeInitialConfigs(QVector<Rectangle*>, QVector <QVector <QPointF> >*) const;
    QVector <Rectangle*> makeConfigC(Rectangle*, QVector <QVector <QPointF> >*) const;
    QVector <Rectangle*> makeConfigR(Rectangle*, QVector <QVector <QPointF> >*) const;
    QVector <Rectangle*> makeConfigT(Rectangle*, QVector <QVector <QPointF> >*) const;
    QVector <Rectangle*> makeConfigRS(Rectangle*, Rectangle*, QVector <QVector <QPointF> >*) const;
    QVector <Rectangle*> makeConfigRA(Rectangle*, Rectangle*, QVector <QVector <QPointF> >*) const;
    void findFeasiblePlacement(Rectangle *, QVector<Rectangle*> *, QVector<QPointF> *) const;
    bool intersectsCircle(Rectangle*) const;
    void processCaseA(QVector <QVector <QPointF>> *, QVector<Rectangle*>*, QVector<QPointF>*, Rectangle*, qreal, qreal, QPointF*, QPointF*, qreal, qreal, QPointF*, QPointF*, QPointF*, int) const;
    void processCaseB(QVector <QVector <QPointF>> *, QVector<Rectangle*>*, QVector<QPointF>*, Rectangle*, qreal, qreal, QPointF*, QPointF*, qreal, qreal, QPointF*, QPointF*, QPointF*, int) const;
    void processCaseC(QVector <QVector <QPointF>> *, QVector<Rectangle*>*, QVector<QPointF>*, Rectangle*, qreal, qreal, QPointF*, QPointF*, qreal, qreal, QPointF*, QPointF*, QPointF*, int) const;
    void processCaseD(QVector <QVector <QPointF>> *, QVector<Rectangle*>*, QVector<QPointF>*, Rectangle*, qreal, qreal, QPointF*, QPointF*, qreal, qreal, QPointF*, QPointF*, QPointF*, int) const;
    void processCaseE(QVector <QVector <QPointF>> *, QVector<Rectangle*>*, QVector<QPointF>*, Rectangle*, qreal, qreal, QPointF*, QPointF*, qreal, qreal, QPointF*, QPointF*, QPointF*, int) const;
    qreal calculateLambdaHorizontal(QPointF *, qreal, QPointF*, qreal, qreal, QPointF*) const;
    qreal calculateLambdaVertical(QPointF *, qreal, QPointF*, qreal, qreal, QPointF*) const;
    int findNearestToCenter(QVector <Rectangle*> *) const;
    qreal calcDistance(QPointF) const;
    qreal calcP(QVector<Rectangle*> *) const;
    QVector<Rectangle*> shiftOnScreen(QVector<Rectangle*>) const;

private:
    QVector<Rectangle*> _rectangles;
    QPointF _circleCenter;
    QPointF _screenCenter;
    QVector <QVector <Rectangle*> > _configurations;
    QVector <Rectangle*> _res;
    QVector <Rectangle*> _finalRes;
    int _poluprecnik;
};
#endif // RECTANGLEPACKING_H
