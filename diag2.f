      SUBROUTINE DIAG2(M,N,D,X)                                         ABHV1480
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ABHV1481
      PARAMETER (MAXDIM=600)                                            ABHV1482
      DIMENSION   D(M), X(M,M)                                          ABHV1483
      DIMENSION   E(MAXDIM)                                             ABHV1484
      IF(M.GT.MAXDIM) THEN                                              ABHV1485
        WRITE(6,10) M,MAXDIM                                            ABHV1486
10      FORMAT(' DIMENSION TOO LARGE IN DIAG2:',2I8)                    ABHV1487
        STOP                                                            ABHV1488
      END IF                                                            ABHV1489
CSTART VAX UNIX-HP UNIX-CONVEX CRAY IBM UNIVAC                          ABHV1490
      EPS=2.5E-13                                                       ABHV1491
      TOL=2.5D-18                                                       ABHV1492
CEND                                                                    ABHV1493
C                                                                       ABHV1494
      IF(N.EQ.1) GO TO 400                                              ABHV1495
C                                                                       ABHV1496
C     HOUSEHOLDER'S REDUCTION                                           ABHV1497
C     SIMULATION OF LOOP DO 150 I=N,2,(-1)                              ABHV1498
C                                                                       ABHV1499
      DO 150 NI=2,N                                                     ABHV1500
      I=N+2-NI                                                          ABHV1501
      L=I-2                                                             ABHV1502
      H=0.0                                                             ABHV1503
      G=X(I,I-1)                                                        ABHV1504
      IF(L) 140,140,20                                                  ABHV1505
   20 DO 30 K=1,L                                                       ABHV1506
   30 H=H+X(I,K)**2                                                     ABHV1507
      S=H+G*G                                                           ABHV1508
      IF(S.GE.TOL) GO TO 50                                             ABHV1509
   40 H=0.0                                                             ABHV1510
      GO TO 140                                                         ABHV1511
   50 IF(H) 140,140,60                                                  ABHV1512
   60 L=L+1                                                             ABHV1513
      F=G                                                               ABHV1514
      G=DSQRT(S)                                                        ABHV1515
      IF(F) 75,75,70                                                    ABHV1516
   70 G=-G                                                              ABHV1517
   75 H=S-F*G                                                           ABHV1518
      X(I,I-1)=F-G                                                      ABHV1519
      F=0.0                                                             ABHV1520
C                                                                       ABHV1521
      DO 110 J=1,L                                                      ABHV1522
      X(J,I)=X(I,J)/H                                                   ABHV1523
      S=0.0                                                             ABHV1524
      DO 80 K=1,J                                                       ABHV1525
   80 S=S+X(J,K)*X(I,K)                                                 ABHV1526
      J1=J+1                                                            ABHV1527
      IF(J1.GT.L) GO TO 100                                             ABHV1528
      DO 90 K=J1,L                                                      ABHV1529
   90 S=S+X(K,J)*X(I,K)                                                 ABHV1530
  100 E(J)=S/H                                                          ABHV1531
  110 F=F+S*X(J,I)                                                      ABHV1532
C                                                                       ABHV1533
      F=F/(H+H)                                                         ABHV1534
C                                                                       ABHV1535
      DO 120 J=1,L                                                      ABHV1536
  120 E(J)=E(J)-F*X(I,J)                                                ABHV1537
C                                                                       ABHV1538
      DO 130 J=1,L                                                      ABHV1539
      F=X(I,J)                                                          ABHV1540
      S=E(J)                                                            ABHV1541
      DO 130 K=1,J                                                      ABHV1542
  130 X(J,K)=X(J,K)-F*E(K)-X(I,K)*S                                     ABHV1543
C                                                                       ABHV1544
  140 D(I)=H                                                            ABHV1545
  150 E(I-1)=G                                                          ABHV1546
C                                                                       ABHV1547
C     ACCUMULATION OF TRANSFORMATION MATRICES                           ABHV1548
C                                                                       ABHV1549
  160 D(1)=X(1,1)                                                       ABHV1550
      X(1,1)=1.0                                                        ABHV1551
      DO 220 I=2,N                                                      ABHV1552
      L=I-1                                                             ABHV1553
      IF(D(I)) 200,200,170                                              ABHV1554
  170 DO 190 J=1,L                                                      ABHV1555
      S=0.0                                                             ABHV1556
      DO 180 K=1,L                                                      ABHV1557
  180 S=S+X(I,K)*X(K,J)                                                 ABHV1558
      DO 190 K=1,L                                                      ABHV1559
  190 X(K,J)=X(K,J)-S*X(K,I)                                            ABHV1560
  200 D(I)=X(I,I)                                                       ABHV1561
      X(I,I)=1.0                                                        ABHV1562
  210 DO 220 J=1,L                                                      ABHV1563
      X(I,J)=0.0                                                        ABHV1564
  220 X(J,I)=0.0                                                        ABHV1565
C                                                                       ABHV1566
C     DIAGONALIZATION OF THE TRIDIAGONAL MATRIX                         ABHV1567
C                                                                       ABHV1568
      B=0.0                                                             ABHV1569
      F=0.0                                                             ABHV1570
      E(N)=0.0                                                          ABHV1571
C                                                                       ABHV1572
      DO 340 L=1,N                                                      ABHV1573
      H=EPS*(DABS(D(L))+DABS(E(L)))                                     ABHV1574
      IF (H.GT.B) B=H                                                   ABHV1575
C                                                                       ABHV1576
C     TEST FOR SPLITTING                                                ABHV1577
C                                                                       ABHV1578
      DO 240 J=L,N                                                      ABHV1579
      IF (DABS(E(J)).LE.B) GOTO 250                                     ABHV1580
  240 CONTINUE                                                          ABHV1581
C                                                                       ABHV1582
C     TEST FOR CONVERGENCE                                              ABHV1583
C                                                                       ABHV1584
  250 IF(J.EQ.L) GO TO 340                                              ABHV1585
C                                                                       ABHV1586
C     SHIFT FROM UPPER 2*2 MINOR                                        ABHV1587
C                                                                       ABHV1588
  260 P=(D(L+1)-D(L))*0.5/E(L)                                          ABHV1589
      R=DSQRT(P*P+1.0)                                                  ABHV1590
      IF(P) 270,280,280                                                 ABHV1591
  270 P=P-R                                                             ABHV1592
      GO TO 290                                                         ABHV1593
  280 P=P+R                                                             ABHV1594
  290 H=D(L)-E(L)/P                                                     ABHV1595
      DO 300 I=L,N                                                      ABHV1596
  300 D(I)=D(I)-H                                                       ABHV1597
      F=F+H                                                             ABHV1598
C                                                                       ABHV1599
C     QR TRANSFORMATION                                                 ABHV1600
C                                                                       ABHV1601
      P=D(J)                                                            ABHV1602
      C=1.0                                                             ABHV1603
      S=0.0                                                             ABHV1604
C                                                                       ABHV1605
C     SIMULATION OF LOOP DO 330 I=J-1,L,(-1)                            ABHV1606
C                                                                       ABHV1607
      J1=J-1                                                            ABHV1608
      DO 330 NI=L,J1                                                    ABHV1609
      I=L+J1-NI                                                         ABHV1610
      G=C*E(I)                                                          ABHV1611
      H=C*P                                                             ABHV1612
C                                                                       ABHV1613
C     PROTECTION AGAINST UNDERFLOW OF EXPONENTS                         ABHV1614
C                                                                       ABHV1615
      IF (DABS(P).LT.DABS(E(I))) GOTO 310                               ABHV1616
      C=E(I)/P                                                          ABHV1617
      R=DSQRT(C*C+1.0)                                                  ABHV1618
      E(I+1)=S*P*R                                                      ABHV1619
      S=C/R                                                             ABHV1620
      C=1.0/R                                                           ABHV1621
      GO TO 320                                                         ABHV1622
  310 C=P/E(I)                                                          ABHV1623
      R=DSQRT(C*C+1.0)                                                  ABHV1624
      E(I+1)=S*E(I)*R                                                   ABHV1625
      S=1.0/R                                                           ABHV1626
      C=C/R                                                             ABHV1627
  320 P=C*D(I)-S*G                                                      ABHV1628
      D(I+1)=H+S*(C*G+S*D(I))                                           ABHV1629
      DO 330 K=1,N                                                      ABHV1630
      H=X(K,I+1)                                                        ABHV1631
      X(K,I+1)=X(K,I)*S+H*C                                             ABHV1632
  330 X(K,I)=X(K,I)*C-H*S                                               ABHV1633
C                                                                       ABHV1634
      E(L)=S*P                                                          ABHV1635
      D(L)=C*P                                                          ABHV1636
      IF (DABS(E(L)).GT.B) GO TO 260                                    ABHV1637
C                                                                       ABHV1638
C     CONVERGENCE                                                       ABHV1639
C                                                                       ABHV1640
  340 D(L)=D(L)+F                                                       ABHV1641
C                                                                       ABHV1642
C     ORDERING OF EIGENVALUES                                           ABHV1643
C                                                                       ABHV1644
      NI=N-1                                                            ABHV1645
  350 DO 380I=1,NI                                                      ABHV1646
      K=I                                                               ABHV1647
      P=D(I)                                                            ABHV1648
      J1=I+1                                                            ABHV1649
      DO 360J=J1,N                                                      ABHV1650
      IF(D(J).GE.P) GOTO 360                                            ABHV1651
      K=J                                                               ABHV1652
      P=D(J)                                                            ABHV1653
  360 CONTINUE                                                          ABHV1654
      IF (K.EQ.I) GOTO 380                                              ABHV1655
      D(K) =D(I)                                                        ABHV1656
      D(I)=P                                                            ABHV1657
      DO 370 J=1,N                                                      ABHV1658
      P=X(J,I)                                                          ABHV1659
      X(J,I)=X(J,K)                                                     ABHV1660
  370 X(J,K)=P                                                          ABHV1661
  380 CONTINUE                                                          ABHV1662
C                                                                       ABHV1663
C     FIXING OF SIGN                                                    ABHV1664
C                                                                       ABHV1665
      DO 385 I=1,N                                                      ABHV1666
      PM=0                                                              ABHV1667
      DO 386 J=1,N                                                      ABHV1668
      IF(PM.GT.DABS(X(J,I))) GOTO 386                                   ABHV1669
      PM =DABS(X(J,I))                                                  ABHV1670
      K=J                                                               ABHV1671
  386 CONTINUE                                                          ABHV1672
      IF(X(K,I).GE.0) GOTO 385                                          ABHV1673
      DO 387 J=1,N                                                      ABHV1674
  387 X(J,I)=-X(J,I)                                                    ABHV1675
  385 CONTINUE                                                          ABHV1676
  390 GO TO 410                                                         ABHV1677
C                                                                       ABHV1678
C     SPECIAL TREATMENT OF CASE N = 1                                   ABHV1679
C                                                                       ABHV1680
  400 D(1)=X(1,1)                                                       ABHV1681
      X(1,1)=1.0                                                        ABHV1682
  410 RETURN                                                            ABHV1683
      END SUBROUTINE                                                    ABHV1684
