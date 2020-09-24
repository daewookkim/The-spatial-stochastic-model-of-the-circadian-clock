breed [ cells cell ]
breed [ nuclei nucleus ]

breed [ obstacles obstacle ]
breed [ mRNArepressors mRNArepressor ]
breed [ repressors repressor ]
breed [ phosrepressors phosrepressor ]

breed [ frepressors frepressor ]
breed [ fphosrepressors fphosrepressor ]


globals [

  counting

  ;; Global variables, which are used to set the size of the cell and that of the nucleus
  cellsize
  nuclradi

  ;; Global variables about the number of activators in the local regions for the nucleus (Fig S9A).
  area1
  area2
  anum1
  anum2
  acon


  ;; Global variables about the reation probabilities for hyperphosphorylation and dephosphorylation
  phosdata
  dephosdata
  phosmatrix
  dephosmatrix

  ;; Global variables, which are used to calcualte the local concentrations of hypophos. PER and hyperphos. PER.
  rad
  ssize
  volume

  ;; Global variable, which is used to determine how the cytoplasmic obstacles decrease the speed of the PER protein.
  ;; E.g. slow=0.5 indicates that PER protein speed halved.
  slow
  ;; Global variable, which describe the speed of Per mRNA, that of PER protein and that of obstacle
  speedp
  speedo

]

extensions [csv matrix]

to load

  clear-all

  set phosdata csv:from-file "phosphorylation_reaction_probability.csv"
  set dephosdata csv:from-file "dephosphorylation_reaction_probability.csv"

  set phosmatrix matrix:from-row-list phosdata
  set dephosmatrix matrix:from-row-list dephosdata
end

to setup

  clear-ticks
  clear-turtles
  clear-patches
  clear-drawing
  clear-all-plots
  clear-output

  draw-box

  set counting 0;

  set cellsize 25
  set nuclradi (cellsize / 6);

  add1 cells 1
  add1 nuclei 1
  add2 obstacles O
  add2 repressors Pc
  add3 mRNArepressors M
  add3 phosrepressors P

  set area1 ((cellsize / 12) ^ 2 * pi) / 8
  set area2 ((cellsize / 6) ^ 2 * pi) / 8 - ((cellsize / 12) ^ 2 * pi) / 8
  set anum1 round (Atot * (area1 / ((cellsize / 6) ^ 2 * pi)))
  set anum2 round (Atot * (area2 / ((cellsize / 6) ^ 2 * pi)))
  set acon (Atot / ((cellsize / 6) ^ 2 * pi))

  set rad 2
  set ssize (rad - 0.25)
  set volume 10

  set slow 0.5

  set speedp (Dp * (cellsize / 2))
  set speedo (Do * (cellsize / 2))

  reset-ticks

end

to go

  ask repressors [

  hatch-frepressors 1 [setshape]

  ]

  ask phosrepressors [

  hatch-fphosrepressors 1 [setshape]

  ]

  ask turtles with [breed = repressors or breed = phosrepressors] [

  ifelse breed = phosrepressors [

    ifelse distancexy 0 0 >= (cellsize / 6) [

    let num (count frepressors with [(distance myself <= ssize) and ((distancexy 0 0) >= (cellsize / 6))])
    let pnum (count fphosrepressors with [(distance myself <= ssize) and ((distancexy 0 0) >= (cellsize / 6))])

    if distancexy 0 0 < (cellsize / 6  + rad) [
     let d distancexy  0 0
     let xstar ((nuclradi ^ 2 - rad ^ 2 + d ^ 2) / (2 * d))
     let ystar sqrt (nuclradi ^ 2 - xstar ^ 2)
     let areanorm (((rad ^ 2) * Pi) - (nuclradi ^ 2 * (asin (ystar / nuclradi) * (pi / 180)) + rad ^ 2 * (asin (ystar / rad) * (pi / 180)) - ystar * (xstar + sqrt (rad ^ 2 - nuclradi ^ 2 + xstar ^ 2))))

     set num (round ((((rad ^ 2) * Pi) / areanorm) * num))
     set pnum (round ((((rad ^ 2) * Pi) / areanorm) * pnum))
    ]

    if pnum > 0 [

     let paseprob ((matrix:get dephosmatrix num pnum) / (pnum / volume))

     if random-float 1 > (e ^ (-(paseprob + (pd3)))) [

        ifelse random-float 1 <= ((pd3) / (paseprob + (pd3)))
        [
        die
        ]
        [
        set breed repressors
        setshape


        ]
       ]
      ]


    ]
    [
      if random-float 1 > (e ^ (-(pd3))) [
         die
       ]
    ]


  ]
  [
    ifelse distancexy 0 0 >= (cellsize / 6)
    [
    let num (count frepressors with [(distance myself <= ssize) and ((distancexy 0 0) >= (cellsize / 6))])
    let pnum (count fphosrepressors with [(distance myself <= ssize) and ((distancexy 0 0) >= (cellsize / 6))])
    let tot (num + pnum)
    let ratio (pnum / tot)

    if distancexy 0 0 < (cellsize / 6  + rad) [
     let d distancexy  0 0
     let xstar ((nuclradi ^ 2 - rad ^ 2 + d ^ 2) / (2 * d))
     let ystar sqrt (nuclradi ^ 2 - xstar ^ 2)
     let areanorm (((rad ^ 2) * Pi) - (nuclradi ^ 2 * (asin (ystar / nuclradi) * (pi / 180)) + rad ^ 2 * (asin (ystar / rad) * (pi / 180)) - ystar * (xstar + sqrt (rad ^ 2 - nuclradi ^ 2 + xstar ^ 2))))
     set num (round ((((rad ^ 2) * Pi) / areanorm) * num))
     set pnum (round ((((rad ^ 2) * Pi) / areanorm) * pnum))

    ]

    if num > 0 [

      let phosprob ((matrix:get phosmatrix num pnum) / (num / volume))

      if random-float 1 > (e ^ (-(phosprob + (pd2)))) [

         ifelse random-float 1 <= ((pd2) / (phosprob + (pd2)))
         [
         die
         ]
         [
         set breed phosrepressors
         setshape

         ]
       ]
      ]
    ]
    [

      if random-float 1 > (e ^ (-(pd2))) [
         die
       ]
    ]

  ]

  ]

ask phosrepressors [
    phosre_active_transport
  ]

ask repressors [
    re_active_transport
  ]

ask mRNArepressors [

    ifelse ( sqrt (xcor * xcor + ycor * ycor) - size / 2 >= (cellsize / 3) )
    [
     if random-float 1 > (e ^ ( - ((pa2 + pd1)))) [

       ifelse random-float 1 <= (pa2 / (pa2 + pd1)) [
         translation
       ]
       [
         die
       ]
      ]
    ]
    [
     if random-float 1 > (e ^ (- (pd1))) [
        die
      ]
    ]

   mRNA_movement

  ]

  transcription11
  transcription12
  transcription13
  transcription14
  transcription15
  transcription16
  transcription17
  transcription18

  transcription21
  transcription22
  transcription23
  transcription24
  transcription25
  transcription26
  transcription27
  transcription28

  ask obstacles [

    obstacle_movement
  ]

  ask frepressors [
   die
  ]

  ask fphosrepressors [
   die
  ]

  tick

  set counting (counting + 1)

   ;;file-open (word "mRNA.txt")
   ;;file-write (count mRNArepressors)
   ;;file-close-all

   file-open (word "totalPER.txt")
   file-write (count repressors) + (count phosrepressors)
   file-close-all

   file-open (word "HypophosPER.txt")
   file-write (count repressors)
   file-close-all

   file-open (word "HyperphosPER.txt")
   file-write (count phosrepressors)
   file-close-all

end

to move2

  ifelse any? obstacles with [overlap myself > 0 ]
  [
    rt random-float 360
    fd speedp * slow

  ]
  [
    rt random-float 360
    fd speedp
  ]

end

to phosre_active_transport

  ifelse sqrt (xcor * xcor + ycor * ycor) > cellsize / 6
  [

  ifelse any? obstacles with [overlap myself > 0]
  [
    let x1 xcor
    let y1 ycor

    rt random-float 360
    fd speedp * slow

    let x2 xcor
    let y2 ycor

    if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
     [
      ifelse random-float 1 > pim
      [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt


        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]
      ]
      [
        let thetha atan x2 y2
        setxy ((cellsize / 6 - size / 2) * (sin thetha)) ((cellsize / 6 - size / 2) * (cos thetha))
      ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt


        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]

  ]
  [
    ifelse (random-float 1 < padvec)
    [

      let x1 xcor
      let y1 ycor

      set heading towardsxy 0 0
      fd speedp

      let x2 xcor
      let y2 ycor

      if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
       [
        ifelse random-float 1 > pim
        [
        let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
        let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
        let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
        let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
        let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

        let xt1 (t1 * x1 + (1 - t1) * x2)
        let yt1 (t1 * y1 + (1 - t1) * y2)

        let xt2 (t2 * x1 + (1 - t2) * x2)
        let yt2 (t2 * y1 + (1 - t2) * y2)


        let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
        let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

        ifelse dis1 < dis2
          [

           let xt (t1 * x1 + (1 - t1) * x2)
           let yt (t1 * y1 + (1 - t1) * y2)

           setxy xt yt

          ]
          [

           let xt (t2 * x1 + (1 - t2) * x2)
           let yt (t2 * y1 + (1 - t2) * y2)
           setxy xt yt

          ]
        ]
        [
          let thetha atan x2 y2
          setxy ((cellsize / 6 - size / 2) * (sin thetha)) ((cellsize / 6 - size / 2) * (cos thetha))
        ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt
        ]

     ]

    ]
    [

      let x1 xcor
      let y1 ycor

      rt random-float 360
      fd speedp

      let x2 xcor
      let y2 ycor

      if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
       [
        ifelse random-float 1 > pim
        [
        let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
        let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
        let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
        let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
        let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

        let xt1 (t1 * x1 + (1 - t1) * x2)
        let yt1 (t1 * y1 + (1 - t1) * y2)

        let xt2 (t2 * x1 + (1 - t2) * x2)
        let yt2 (t2 * y1 + (1 - t2) * y2)


        let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
        let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

        ifelse dis1 < dis2
          [
           let xt (t1 * x1 + (1 - t1) * x2)
           let yt (t1 * y1 + (1 - t1) * y2)
           setxy xt yt

          ]
          [
           let xt (t2 * x1 + (1 - t2) * x2)
           let yt (t2 * y1 + (1 - t2) * y2)
           setxy xt yt

          ]
        ]
        [
          let thetha atan x2 y2
          setxy ((cellsize / 6 - size / 2) * (sin thetha)) ((cellsize / 6 - size / 2) * (cos thetha))
        ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt
        ]

     ]


    ]
  ]
  ]
  [

    let x1 xcor
    let y1 ycor

    rt random-float 360
    fd speedp

    let x2 xcor
    let y2 ycor

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 6
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]

  ]

end

to re_active_transport

  ifelse any? obstacles with [overlap myself > 0]
  [
    let x1 xcor
    let y1 ycor

    rt random-float 360
    fd speedp * slow

    let x2 xcor
    let y2 ycor

    if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

        ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt
        ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]

  ]
  [
    ifelse (random-float 1 < padvec)
    [

      let x1 xcor
      let y1 ycor

      set heading towardsxy 0 0
      fd speedp

      let x2 xcor
      let y2 ycor

      if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
       [
        let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
        let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
        let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
        let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
        let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

        let xt1 (t1 * x1 + (1 - t1) * x2)
        let yt1 (t1 * y1 + (1 - t1) * y2)

        let xt2 (t2 * x1 + (1 - t2) * x2)
        let yt2 (t2 * y1 + (1 - t2) * y2)

        let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
        let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

        ifelse dis1 < dis2
          [
           let xt (t1 * x1 + (1 - t1) * x2)
           let yt (t1 * y1 + (1 - t1) * y2)
           setxy xt yt

          ]
          [
           let xt (t2 * x1 + (1 - t2) * x2)
           let yt (t2 * y1 + (1 - t2) * y2)
           setxy xt yt
          ]
     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]

    ]
    [

      let x1 xcor
      let y1 ycor

      rt random-float 360
      fd speedp

      let x2 xcor
      let y2 ycor

      if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
       [
        let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
        let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
        let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
        let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
        let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

        let xt1 (t1 * x1 + (1 - t1) * x2)
        let yt1 (t1 * y1 + (1 - t1) * y2)

        let xt2 (t2 * x1 + (1 - t2) * x2)
        let yt2 (t2 * y1 + (1 - t2) * y2)

        let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
        let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

        ifelse dis1 < dis2
          [
           let xt (t1 * x1 + (1 - t1) * x2)
           let yt (t1 * y1 + (1 - t1) * y2)
           setxy xt yt

          ]
          [
           let xt (t2 * x1 + (1 - t2) * x2)
           let yt (t2 * y1 + (1 - t2) * y2)
           setxy xt yt

          ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt

        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]



    ]
  ]
end

to mRNA_movement

  ifelse sqrt (xcor * xcor + ycor * ycor) > cellsize / 6
  [
    let x1 xcor
    let y1 ycor

    rt random-float 360
    fd speedp

    let x2 xcor
    let y2 ycor

    if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
     [

      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt


        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt
        ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [

      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt


        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]


  ]
  [
    rt random-float 360
    fd speedp
  ]

end


to obstacle_movement

    let x1 xcor
    let y1 ycor

    rt random-float 360
    fd speedo

    let x2 xcor
    let y2 ycor

    if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
     [

      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt


        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt
        ]

     ]

    if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
     [

      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
         let xt (t1 * x1 + (1 - t1) * x2)
         let yt (t1 * y1 + (1 - t1) * y2)
         setxy xt yt


        ]
        [
         let xt (t2 * x1 + (1 - t2) * x2)
         let yt (t2 * y1 + (1 - t2) * y2)
         setxy xt yt

        ]

     ]

end

to-report overlap [ agent ]
   report (size + [size] of agent) / 2 - distance agent
end

;; reactions

to translation

  hatch-repressors 1 [setshape]

end

to transcription11
  let recon1 (count fphosrepressors with [ (distancexy 0 0) = 0]) / area1
  let recon2 (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (0 <= (atan xcor ycor)) and ((atan xcor ycor) < 360 / 8)]) / area1
  let recon (recon1 + recon2)
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt11 mRNArepressors transnum

end

to transcription12

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 2 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt12 mRNArepressors transnum


end

to transcription13

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (2 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 3 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt13 mRNArepressors transnum

end

to transcription14

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (3 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 4 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt14 mRNArepressors transnum

end

to transcription15

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (4 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 5 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt15 mRNArepressors transnum

end

to transcription16

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (5 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 6 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt16 mRNArepressors transnum

end

to transcription17

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (6 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 7 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt17 mRNArepressors transnum

end

to transcription18

  let recon (count fphosrepressors with [((distancexy 0 0) > 0) and ((distancexy 0 0) < cellsize / 12) and (7 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 8 * 360 / 8)]) / area1
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum1 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt18 mRNArepressors transnum

end

to transcription21

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (0 <= (atan xcor ycor)) and ((atan xcor ycor) < 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt21 mRNArepressors transnum

end

to transcription22

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 2 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt22 mRNArepressors transnum

end

to transcription23

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (2 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 3 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt23 mRNArepressors transnum


end

to transcription24

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (3 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 4 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt24 mRNArepressors transnum


end

to transcription25

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (4 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 5 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt25 mRNArepressors transnum

end

to transcription26

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (5 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 6 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt26 mRNArepressors transnum

end

to transcription27

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (6 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 7 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt27 mRNArepressors transnum

end

to transcription28

  let recon (count fphosrepressors with [(cellsize / 12 <= (distancexy 0 0)) and ((distancexy 0 0) < cellsize / 6) and (7 * 360 / 8 <= (atan xcor ycor)) and ((atan xcor ycor) < 8 * 360 / 8)]) / area2
  let transf ((acon - recon - Kd + sqrt ( (acon - recon - Kd) ^ 2 + 4 * acon * Kd )) / (2 * acon))
  let acounting 0
  let transnum 0
  while [ acounting <= (round (anum2 * transf)) ] [

    if random-float 1 > (e ^ (- (pa1))) [

      set transnum (transnum + 1)

    ]
  set acounting (acounting + 1)
  ]

  addt28 mRNArepressors transnum

end

to addt11 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha random-float 360 / 8
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt12 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((360 / 8) + (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt13 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((2 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt14 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((3 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt15 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((4 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt16 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((5 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt17 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((6 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt18 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 12)
      let thetha ((7 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt21 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha random-float 360 / 8
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt22 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((360 / 8) + (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt23 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((2 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt24 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((3 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt25 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((4 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt26 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((5 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt27 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((6 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

to addt28 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius ((cellsize / 12) + random-float (cellsize / 12))
      let thetha ((7 * 360 / 8) +  (random-float 360 / 8))
      setxy (radius * (sin thetha)) (radius * (cos thetha))]
end

;; construct the environments

to draw-box

ask patches [set pcolor white]

end

to add1 [kind amount]
  create-turtles amount
  [ set breed kind
    setshape
    setxy 0 0]
end

to add2 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape

      let radius (cellsize / 6) + size / 2 + random-float (cellsize / 3 - size)
      let thetha random-float 360
      setxy (radius * (cos thetha)) (radius * (sin thetha))]

end

to add3 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 2) - size / 2
      let thetha random-float 360

      let xvalue (radius * (cos thetha))
      let yvalue (radius * (sin thetha))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)

      while [dist + size / 2 > cellsize / 6 and dist - size / 2 < cellsize / 6]
      [
      set radius random-float (cellsize / 2) - size / 2
      set thetha random-float 360

      set xvalue (radius * (cos thetha))
      set yvalue (radius * (sin thetha))
      set dist sqrt (xvalue * xvalue + yvalue * yvalue)

      ]
      setxy xvalue yvalue]
end

;; procedure that assigns a specific shape to a turtle

to setshape

  if breed = cells
 [ set color White
   set shape "circle"
   set size cellsize  ]

 if breed = nuclei
 [set color grey
  set shape "circle"
  set color lput 10 extract-rgb color
  set size (cellsize / 3)  ]

 if breed = obstacles
 [set color black
  set color lput 100 extract-rgb color
  set shape "circle"
  set size 1]

  if breed = mRNArepressors
 [set color cyan
  set shape "circle"
  set size 0.5]

  if breed = repressors
 [set color Orange
  set shape "circle"
  set size 0.5]

  if breed = phosrepressors
 [set color Violet
  set shape "circle"
  set size 0.5]

  if breed = frepressors
 [set color blue
  set color lput 0 extract-rgb color
  set shape "circle"
  set size 0.5]

  if breed = fphosrepressors
 [set color green
  set color lput 0 extract-rgb color
  set shape "circle"
  set size 0.5]

end
@#$#@#$#@
GRAPHICS-WINDOW
560
85
830
356
-1
-1
8.452
1
10
1
1
1
0
1
1
1
-15
15
-15
15
1
1
1
ticks
30.0

BUTTON
153
11
259
53
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
270
10
376
52
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
15
397
862
600
Concentration
time
Conc.
0.0
50000.0
0.0
1300.0
true
true
"" ""
PENS
"mRNA" 1.0 0 -8990512 true "" "plot count mRNArepressors"
"Total" 1.0 0 -16777216 true "" "plot ((count repressors) + (count phosrepressors))"
"Hypo" 1.0 0 -817084 true "" "plot count repressors"
"Hyper" 1.0 0 -6917194 true "" "plot count phosrepressors"

TEXTBOX
15
63
165
81
Initial condition\n
14
0.0
1

SLIDER
7
243
179
276
O
O
0
600
150.0
1
1
NIL
HORIZONTAL

SLIDER
10
88
182
121
Atot
Atot
0
600
500.0
1
1
NIL
HORIZONTAL

SLIDER
9
128
181
161
M
M
0
500
303.0
1
1
NIL
HORIZONTAL

SLIDER
8
203
180
236
P
P
0
1000
701.0
1
1
NIL
HORIZONTAL

SLIDER
9
166
181
199
Pc
Pc
0
1000
555.0
5
1
NIL
HORIZONTAL

TEXTBOX
195
63
270
81
Parameter
14
0.0
1

SLIDER
194
89
366
122
pa1
pa1
0
1
0.015
0.001
1
NIL
HORIZONTAL

SLIDER
194
127
366
160
pa2
pa2
0
1
0.005
0.001
1
NIL
HORIZONTAL

SLIDER
192
165
364
198
pd1
pd1
0
1
0.001
0.001
1
NIL
HORIZONTAL

SLIDER
191
202
363
235
pd2
pd2
0
1
0.001
0.001
1
NIL
HORIZONTAL

SLIDER
373
89
545
122
Kd
Kd
0
1
1.0E-7
0.0000001
1
NIL
HORIZONTAL

SLIDER
373
127
545
160
Dp
Dp
0
1
0.032
0.001
1
NIL
HORIZONTAL

SLIDER
372
169
544
202
Do
Do
0
1
0.008
0.001
1
NIL
HORIZONTAL

SLIDER
191
239
363
272
pd3
pd3
0
1
0.002
0.001
1
NIL
HORIZONTAL

SLIDER
373
207
545
240
padvec
padvec
0
1
0.2
0.1
1
NIL
HORIZONTAL

SLIDER
373
244
545
277
pim
pim
0
1
0.01
0.001
1
NIL
HORIZONTAL

BUTTON
25
11
128
54
load
load
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

This model was extended from the previous mathematical model of the mammalian circadian clock (Kim et al., 2012, MSB). To investigate the effect of the crowdedness of cytoplasmic contents on the circadina clock, this model was developed.

The model simulation suggests that the spatial regulation of PER protein, the repressor of the circadian clock, leads to its sharp switch-like phosphorylation and thus sharp nuclear translocation. This nonlinear nuclear entry plays the critial role in generating the robust circadian rhythms. Please see Beesley et al., for details description of the model.


## HOW TO USE IT

Step 1. Please set the initial number of agents.

Atot: Total number of activators.

M: The initial number of Per mRNA.

HYPOPER: The initial number of hypophosphorylated PER.

HYPERPER: The initial number of hyperphosphorylated PER.

O: Number of cytopalsmic obstacles. Please choose the value of 0 to regulate the crowdedness of the cytoplasmic contents.
Examples
- The model with 150 cytoplasmic obstacles is the normal cell.
- The model with 275 cytoplasmic obstacles is the overcrowded cell.
- The model with more than 300 cytoplasmoc obstacles is the extremely overcrowded cell (i.e. adipocyte). Thus, the model do not simulate the rhythmic PER expression.



Step 2. Please set the parameters, which are described below. The current parameter setting makes the model to simulate the rhythmic PER expression of the normal cell.

pa1: Reaction probability for Per mRNA production for each time step (i.e. tick).
pa2: Reaction probability for PER protein translation for each time step.
pd1: Reaction probability for Per mRNA degradation for each time step.
pd2: Reaction probability for hypophos. PER degradtion for each time step. 
pd3: Reaction probability for hyperphos. PER degradtion for each time step.

Kd: Dissociation constant between hyperphos. PER and activator.
D: Movement step size of PermRNA and PER protein for each time step.
Dobs: Movement step size of cytoplasmic obstacles for each time step.
padvec: Probability that the PER protein is advected to the peri-nucleus by the cytoplasmic flux.
pim: Probability that hyperphos.PER in the cytoplasm is imported to the nucleus for each time step.



Step 3. Load the reation probabilities for hyperphosphorylation and dephosphrylation by pushing "load" buttom.



Step 4. Set the initial condition of the model by pushing "setup" buttom.



Step 5. Run the simulation by pushing "go" buttom. 





@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

complex
true
0
Polygon -2674135 true false 76 47 197 150 76 254 257 255 257 47
Polygon -10899396 true false 79 46 198 148 78 254

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

enzyme
true
0
Polygon -2674135 true false 76 47 197 150 76 254 257 255 257 47

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

inhib complex
true
0
Polygon -2674135 true false 76 47 197 150 76 254 257 255 257 47
Polygon -1184463 true false 77 48 198 151 78 253 0 253 0 46

inhibitor
true
0
Polygon -1184463 true false 197 151 60 45 1 45 1 255 60 255

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

substrate
true
5
Polygon -10899396 true true 76 47 197 151 75 256

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
