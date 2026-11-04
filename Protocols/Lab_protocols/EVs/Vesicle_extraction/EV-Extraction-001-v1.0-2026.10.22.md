# Quick Overview

  -----------------------------------------------------------------------
  **Parameter**          **Information**
  ---------------------- ------------------------------------------------
  Protocol ID            EV-Extraction-001

  Version                1.0

  Method Owner           Haining Huang

  Estimated Hands-on     \~4--6 h (excluding the unattended
  Time                   centrifugation and the 2 h ultracentrifugation
                         run)

  Total Duration         \~6.5--9 h of bench time; usually spread over
                         1--3 days with optional overnight 4°C holds

  Difficulty Level       ☐ Basic ☐ Intermediate ☑ Advanced

  Training Required      ☑ Yes ☐ No

  Biosafety Level        Depends on the strain used; in all cases, work
                         must comply with the applicable institutional
                         laboratory safety regulations

  Typical Throughput     Max \~5 L of culture processed each time

  Sample types           Bacterial culture supernatant (e.g., *E. coli*
                         or *P. putida* liquid culture)

  Output                 Purified bacterial extracellular vesicles (EVs),
                         resuspended in \~1 mL PBS

  Related Protocols      EV-Purification-001 (Density Gradient
                         Purification --- see Optional Step)
  -----------------------------------------------------------------------

# Method Information

  ---------------------- ------------------------------------------------
  Method Title           SOP for Bacterial Extracellular Vesicles (EVs)
                         Extraction

  Author(s)              Haining Huang

  Date Created           2026-08-28

  Last Updated           2026-10-21

  Reviewed By            Chika Ejikeugwu

  Approved By            
  ---------------------- ------------------------------------------------

# Purpose

To isolate and enrich extracellular vesicles (EVs) from bacterial
culture supernatant through sequential centrifugation, filtration,
concentration, and ultracentrifugation, yielding an EV preparation
suitable for downstream functional validation assays. If proteomic
analysis is required, further purification by density gradient can be
performed (Optional Step).

### Typical Starting Material

Typically, *E. coli* MG1655 carrying plasmid pKJK5 or *P. putida* KT2440
carrying plasmid RSF1010 is cultured in LB medium at 37°C for 16 h at
150 rpm until it reaches an OD600 of approximately 1.5. Approximately 2
L of culture is processed per batch.

# Scope

### Suitable for

- EV isolation from bacterial liquid culture supernatant (broth culture)

- Sample volumes from \<1 L (dead-end filtration can be used) up to \~5
  L (using the TFF system)

- EV preparations intended for downstream functional validation assays

### Not suitable for

- Non-liquid/solid samples

# Safety Considerations

### Hazards

**Biological:** Bacterial culture supernatant may contain live bacterial
cells or residual pathogens prior to filtration; handle at the
appropriate biosafety level

**Chemical:** 0.01 M NaOH (mildly alkaline, used for TFF system storage)
--- avoid contact with skin and eyes; wear safety goggles

**Flammable:** Ethanol (70--75%, lab-grade; used to clean
ultracentrifuge tubes and filter inserts) --- keep away from open flames
and other ignition sources

**Physical:** Ultracentrifugation at 100,000 × g --- if a tube is
underfilled, the vacuum drawn during the run may crush the thin-wall
tube (see Critical Control Points below); observe standard mechanical
safety when operating the centrifuge/ultracentrifuge

### Required PPE

- ☑ Lab coat

- ☑ Gloves

- ☑ Safety goggles (when handling NaOH or ethanol)

### Waste Disposal

**Biological waste:** Centrifuge bottles and ultracentrifuge tubes,
filter membranes and TFF cassettes that have been in contact with
bacterial culture are decontaminated before reuse or disposal (see
"Decontamination of Reusable Bottles, Tubes and Glassware" above).
Single-use items --- filter membranes, syringe filters and pipette tips
--- are discarded as biological waste; reusable items are collected
centrally for decontamination and washing.

**Chemical waste:** 0.01 M NaOH solution --- dispose of per
institutional chemical waste guidelines

### Decontamination of Reusable Bottles, Tubes and Glassware

Bottles, tubes and glassware that have been in contact with bacterial
culture are decontaminated by **one** of the two routes below before
washing. Use one route only: do not autoclave items that have already
been chemically disinfected, and avoid autoclaving the same item
repeatedly --- the plastic degrades with every cycle.

**Route A --- Autoclaving.** Separate the bottles from their caps, bag
them and autoclave, then wash. This is the most thorough option, but
polypropylene bottles and thin-wall ultracentrifuge tubes tolerate only
a limited number of cycles; keep a record of the cycles for items that
are reused often.

**Route B --- Chemical disinfection.** Immerse the items in 2% (v/v)
Wofasteril (peracetic acid) for 1--2 h, then rinse thoroughly with water
and wash as usual. Use this route for items that would otherwise have to
be autoclaved repeatedly.

# Materials and Reagents

### Equipment

  ------------------------------------------------------------------------
  **Equipment**           **Model**               **Location**
  ----------------------- ----------------------- ------------------------
  1 L centrifuge bottles  Thermo Scientific       E27/K21 in building 4.1
                                                  (UFZ)

  SORVALL LYNX 6000       Thermo Scientific Inc., K21 in building 4.1
  Centrifuge              Waltham, MA, USA        (UFZ)

  Vacuum pump             Laboratory diaphragm    E29 in building 4.1
                          vacuum pump (any model  (UFZ)
                          suitable for            
                          filtration)             

  Vacuum filtration flask e.g., Büchner/side-arm  E29 in building 4.1
                          filter flask, used with (UFZ)
                          the vacuum pump         

  Concentrator tubes      Amicon® Ultra-15, 100   E29 in building 4.1
                          kDa MWCO                (UFZ)

  Tangential Flow         Interchangeable         One in K29 in building
  Filtration (TFF) system cassettes: Sartorius    4.1 (UFZ), one in BTU
                          0.2 µm cassette (Step   
                          2, filtration) and      
                          Millipore 100 kDa       
                          cassette (Step 3,       
                          Method B)               

  Analytical balance      precision 0.01 g        E29 in building 4.1
                                                  (UFZ)

  Ultracentrifuge + rotor Beckman Coulter Optima  Lab 221 in building 4.1
                          XP-80 ultracentrifuge   (UFZ)
                          with SW 32 Ti           
                          swinging-bucket rotor   
                          and buckets (#1--#6)    

  Ultracentrifuge tubes   Open-Top Thinwall PP,   E29 in building 4.1
                          38.5 mL, 25×89 mm, Cat. (UFZ)
                          No. 326823              
  ------------------------------------------------------------------------

### Materials and Reagents

+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| **Reagent/Material** | **Specification**       | **Supplier**                                                                                                                                                           | **Storage**   |
+======================+=========================+========================================================================================================================================================================+===============+
| 1× PBS (stock)       | Not filtered;           | Prepared in-house                                                                                                                                                      | Room          |
|                      | general-purpose stock   |                                                                                                                                                                        | temperature   |
|                      | solution                |                                                                                                                                                                        | (RT)          |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| 1× PBS (filtered,    | Sterile, 0.2 µm         | Prepared in-house                                                                                                                                                      | 4°C           |
| cold)                | filtered, pre-chilled   |                                                                                                                                                                        |               |
|                      | --- used for            |                                                                                                                                                                        |               |
|                      | resuspension steps      |                                                                                                                                                                        |               |
|                      | (e.g., Step 4.7)        |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| 0.01 M NaOH (working | 0.01 M (in DI water)    | Prepared in-house                                                                                                                                                      | RT            |
| solution)            |                         |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| PBS salts            | NaCl, KCl, Na₂HPO₄      | Any supplier                                                                                                                                                           | RT            |
|                      | (anhydrous), KH₂PO₄ --- |                                                                                                                                                                        |               |
|                      | analytical grade; used  |                                                                                                                                                                        |               |
|                      | to prepare the 1× PBS   |                                                                                                                                                                        |               |
|                      | stock (Table 2)         |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| NaOH (solid)         | Pellets, analytical     | Any supplier                                                                                                                                                           | RT            |
|                      | grade; used to prepare  |                                                                                                                                                                        |               |
|                      | the 0.01 M NaOH         |                                                                                                                                                                        |               |
|                      | solution                |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Ethanol              | 70--75% (v/v),          | Any supplier                                                                                                                                                           | RT ---        |
|                      | lab-grade; used to      |                                                                                                                                                                        | flammable;    |
|                      | clean ultracentrifuge   |                                                                                                                                                                        | store per     |
|                      | tubes and concentrator  |                                                                                                                                                                        | institutional |
|                      | filter inserts          |                                                                                                                                                                        | rules         |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Deionized (DI) water | Large volume            | Available in the lab                                                                                                                                                   | RT            |
|                      | consumable; used for    |                                                                                                                                                                        |               |
|                      | rinsing the TFF system, |                                                                                                                                                                        |               |
|                      | concentrator tubes, and |                                                                                                                                                                        |               |
|                      | ultracentrifuge tubes   |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Media bottle         | 5 L, glass              | Any supplier                                                                                                                                                           | RT            |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Erlenmeyer flask     | 1 L, 2 L or 3 L         | Any supplier                                                                                                                                                           | RT            |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Centrifuge tubes     | 50 mL                   | Any supplier                                                                                                                                                           | RT            |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Filter membranes     | PES, 0.8 / 0.45 / 0.2   | Cytiva (Pall) Supor or GVS                                                                                                                                             | RT            |
| (dead-end            | µm, 47 mm diameter ---  |                                                                                                                                                                        |               |
| filtration)          | see Ordering            | (or any supplier)                                                                                                                                                      |               |
|                      | Information             |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Syringe filters      | 0.2 µm, sterile         | [[https://www.sigmaaldrich.com/DE/de/product/mm/slgpr33rs]{.underline}](https://www.sigmaaldrich.com/DE/de/product/mm/slgpr33rs)                                       | RT            |
|                      |                         |                                                                                                                                                                        |               |
|                      |                         | [[https://www.fishersci.de/shop/products/sterile-pes-syringe-filter/15206869]{.underline}](https://www.fishersci.de/shop/products/sterile-pes-syringe-filter/15206869) |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Syringes             | 20 mL / 50 mL           | Any supplier                                                                                                                                                           | RT            |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| EP tubes             | 1.5 mL/2 mL             | Thermo Fisher                                                                                                                                                          | RT            |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+
| General labware      | Pipettes, motorized     | Any supplier                                                                                                                                                           | RT            |
|                      | pipette filler, tube    |                                                                                                                                                                        |               |
|                      | racks, micropipettes    |                                                                                                                                                                        |               |
|                      | (various volumes),      |                                                                                                                                                                        |               |
|                      | sterile syringes        |                                                                                                                                                                        |               |
+----------------------+-------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------+

# Reagent Preparation

### Solution/Buffer 1 --- 1× PBS (RT stock)

**Table 1. Solution/Buffer 1 --- 1× PBS (RT stock)**

  ---------------- ------------------------------------------------------
  Composition      1× PBS, autoclaved (filtration optional)

  Preparation      See Table 2

  Storage          Room temperature (RT)

  Stability        Stable for ≥1 year at RT
  ---------------- ------------------------------------------------------

**Table 2. Composition of 1× PBS (0.01 M, pH 7.2)**

  -----------------------------------------------------------------------
  **Component**                **Amount per 1 L**
  ---------------------------- ------------------------------------------
  NaCl                         8 g

  KCl                          0.2 g

  Na~2~HPO~4~                  1.44 g

  KH~2~PO~4~                   0.24 g

  Solvent                      Double-distilled water (ddH~2~O), make up
                               to 1 L

  pH                           Adjust to 7.2
  -----------------------------------------------------------------------

### Solution/Buffer 2 --- 1× PBS (filtered, cold)

**Table 3. Solution/Buffer 2 --- 1× PBS (filtered, cold)**

  ---------------- ------------------------------------------------------
  Composition      1× PBS, sterile, 0.2 µm filtered

  Preparation      Filter the PBS stock through a 0.2 µm filter

  Storage          4°C

  Stability        Stable
  ---------------- ------------------------------------------------------

### Solution/Buffer 3 --- 0.01 M NaOH

**Table 4. Solution/Buffer 3 --- 0.01 M NaOH**

  ---------------- ------------------------------------------------------
  Composition      0.01 M NaOH in DI water

  Preparation      Dissolve 0.4 g of NaOH in 1 L of DI water.

  Storage          Room temperature (RT)

  Stability        Used to flush/store the TFF system between runs to
                   inhibit bacterial growth; very stable
  ---------------- ------------------------------------------------------

# Procedure

![Workflow overview of the EV extraction
procedure](./images/media/image1.png){width="6.6in" height="6.765in"}

**Figure 1. Workflow overview of the bacterial EV extraction procedure**

*Approximate time for each step is shown on the right. Steps 1--2 remove
bacteria and debris; Steps 3--4 enrich the vesicles.*

### Before Starting --- Checklist

- [ ] Required training completed

- [ ] Reagents prepared (1× PBS, 0.01 M NaOH)

- [ ] Equipment available (TFF system, ultracentrifuge and rotor
  booked/confirmed)

- [ ] Calibration verified (analytical balance)

- [ ] Sample IDs assigned

### Step 1: Centrifugation (Bacterial Removal)

*Note: Because this is a bacterial culture (unlike, e.g., soil or
environmental samples), the starting density of the culture is
relatively uniform, so balancing between bottles is straightforward.
Bottles must still be balanced carefully before every run (see Critical
Control Points).*

1.  Pour the bacterial culture into 1 L centrifuge bottles.

2.  Strictly balance the bottles using a scale or by matching the liquid
    volume meniscus.

3.  Centrifuge at 5000--10,000 × g (usually 8000 × g), 30--60 min,
    22--25°C (room temperature).

4.  After centrifugation, slowly transfer the supernatant into a clean
    2/3/5 L media bottle.

5.  Separate the centrifuge bottles and caps and decontaminate them by
    Route A or Route B (see "Decontamination of Reusable Bottles, Tubes
    and Glassware"), then wash.

**Expected result:** Clear supernatant free of bacterial pellet, ready
for filtration.

**⏱ This step will take approximately 1.5 h (including loading,
balancing, centrifugation, and supernatant transfer). If you cannot
finish on the same day, the supernatant can be kept at 4°C and the
protocol continued the next day.**

### Step 2: Filtration (Removal of Residual Bacteria and Debris)

*Note: Traditional dead-end filtration (sequentially using 0.8, 0.45,
and 0.2 µm filters) can be used; for larger volumes, a tangential flow
filtration (TFF) system is recommended to accelerate processing.
However, not all labs have the 0.2 µm cassette.*

**Sequential Dead-End Filtration**

1.  Assemble the vacuum filtration unit: fit a 47 mm 0.8 µm PES membrane
    into the filter holder and connect the filtration flask to the
    vacuum pump. The membranes do not need to be pre-wetted.

2.  Pour the supernatant from Step 1 into the funnel and filter it under
    vacuum. Replace the membrane whenever the flow rate drops
    noticeably.

3.  Repeat the filtration with a fresh 47 mm 0.45 µm PES membrane, then
    with a fresh 47 mm 0.2 µm PES membrane, collecting the filtrate each
    time. Use a new membrane at every stage.

4.  Discard the used membranes as biological waste and send the
    filtration glassware for sterilization and washing.

**Tangential Flow Filtration (TFF) --- Sartorius 0.2 µm cassette**

1.  Rinse the internal system of the TFF equipment with 2 L of DI water
    to remove residual liquid.

2.  Insert the feed tube into the sample bottle. Collect the permeate;
    route the retentate back into the feed tube to circulate until all
    liquid is collected as permeate.

3.  Flush the lines with 2 L of DI water, then circulate 1 L of 0.01 M
    NaOH to inhibit bacterial growth for storage.

4.  Decontaminate the original sample bottles by Route A or Route B (see
    "Decontamination of Reusable Bottles, Tubes and Glassware"), then
    wash.

**Expected result:** Filtrate free of residual bacteria and debris
(clarified through 0.2 µm).

**⏱ This step will take approximately 1--2 h, depending on sample volume
and filter flux. If you cannot finish on the same day, the filtrate can
be kept at 4°C and the protocol continued the next day.**

**Purpose of Steps 1--2: to remove bacteria and large particles/debris
from the sample.**

### Step 3: Concentration (Volume Reduction)

*Two methods can be used for this step. Both give an equivalent result
--- use whichever system is available:*

**Method A --- Tube-based centrifugal concentration (Amicon® Ultra-15,
100 kDa MWCO)**

*Performed in E29 (building 4.1, UFZ) using an Eppendorf centrifuge and
100 kDa Millipore (Amicon® Ultra-15) concentrator tubes.*

1.  Rinse the 15 mL concentrator tubes with DI water and place the
    filter insert into the plastic tube.

2.  Use a pipette to transfer the filtered sample into the concentrator
    tubes (max capacity per run: 16 tubes × 15 mL = 240 mL).

3.  Centrifuge at 4°C, 2000 × g, 3--5 min (adjust the time based on the
    residual volume in individual tubes).

4.  When the retentate volume is \<1 mL, transfer it to a 50 mL
    centrifuge tube for pooling. Repeat this process until all samples
    are concentrated (est. time: \~2 hours for 2 L).

5.  After the last collection, rinse the concentrator tubes with 1--20
    mL of PBS to recover vesicles attached to the membrane and add the
    rinse to the pooled sample. The pooled volume does not have to fill
    the ultracentrifuge tubes --- they are topped up with PBS in Step
    4.2.

6.  Store the concentrated samples at 4°C.

7.  Rinse the filter inserts multiple times with DI water and store them
    in 0.01 M NaOH or 75% ethanol; wash and dry the plastic outer tubes.

**Expected result:** Retentate volume suitable for loading into the
ultracentrifuge tubes in the next step.

**⏱ Method A will take approximately 2 h for 2 L of sample. If you
cannot finish on the same day, the concentrated sample can be kept at
4°C and the protocol continued the next day.**

**Method B --- Tangential Flow Filtration (TFF) concentration (Millipore
100 kDa TFF system)**

*Use this method when a dedicated 100 kDa TFF concentrator system is
available; it is faster than Method A for large volumes.*

1.  Rinse the TFF concentration system with DI water.

2.  Load the filtered sample (from Step 2) into the TFF system
    reservoir.

3.  Circulate the sample through the 100 kDa cassette, discarding the
    permeate, until \~150 mL of concentrated retentate remains.

4.  Rinse the remaining sample out of the system with PBS and collect a
    further \~80 mL, giving \~230 mL of pooled retentate in total.

5.  Flush the system with DI water, then circulate 0.01 M NaOH for
    storage between uses (as in Step 2).

6.  Store the pooled retentate at 4°C until ultracentrifugation.

**Expected result:** Same as Method A, achieved more quickly.

**⏱ Method B will take approximately 1 h for 2 L of sample (faster than
Method A). If you cannot finish on the same day, the concentrated sample
can be kept at 4°C and the protocol continued the next day.**

*Typical final volume prior to ultracentrifugation (either method): up
to 6 × 38.5 mL (38.5 mL = capacity of each ultracentrifuge tube; up to 6
tubes/samples can be processed per ultracentrifugation run).*

### Step 4: Ultracentrifugation (EV Enrichment)

1.  Filter the concentrated sample using a 0.2 µm syringe filter and
    transfer it into the 38.5 mL ultracentrifuge tubes.

2.  Top up each tube with pre-cooled PBS to \~37--38 mL, close to the
    38.5 mL capacity. Underfilled tubes collapse under vacuum --- see
    Critical Control Points.

3.  Balance the ultracentrifuge tubes precisely to 0.01 g using an
    analytical balance.

4.  Place the tubes symmetrically into the SW 32 Ti rotor buckets (e.g.,
    #1 opposite #4, #2 opposite #5). All six buckets must be attached to
    the rotor for every run, even when fewer tubes are loaded.

5.  Set parameters: Rotor/Tube: SW 32 Ti / 38.5 mL; Speed: 100,000 × g;
    Time: 2 h; Temperature: 4°C; Accel: Max; Decel: 7/8 (or set to no
    coast).

6.  Operation monitoring and run: click "Start vacuum" to evacuate the
    chamber; the run begins after 2--3 min and the rotor reaches the set
    speed in \~5 min. After the run, the vacuum is released once the
    speed drops below 2000 rpm; then click "Open vacuum" before removing
    the tubes.

7.  Remove the tubes and pour off the supernatant completely in one
    smooth motion. Add 100--200 µL of pre-cooled PBS to the tube and
    resuspend the pellet by pipetting up and down repeatedly aimed at
    the bottom (final total volume should be \~1 mL). If more than one
    tube was run, pool the resuspended pellets from all tubes into a
    single collection tube and adjust the combined volume to \~1 mL with
    pre-cooled PBS.

8.  Clean the ultracentrifuge tubes for reuse: decontaminate them by
    Route A or Route B (see "Decontamination of Reusable Bottles, Tubes
    and Glassware"), then rinse with DI water, dry and store at room
    temperature. Thin-wall polypropylene tubes tolerate only a limited
    number of autoclave cycles, so Route B is preferred for tubes that
    are reused often.

**Expected result:** EV pellet resuspended in \~1 mL PBS, ready for
downstream storage or assays.

**⏱ This step will take approximately 3--3.5 h (including tube loading,
balancing, the 2 h run, and pellet resuspension). If you cannot start
the ultracentrifugation run on the same day, the pre-filtered
concentrated sample can be kept at 4°C and the run continued the next
day.**

**Purpose of Steps 3--4: to enrich the vesicles.**

In theory, EVs are obtained once the steps above are complete. Depending
on your experimental design, however, consider whether further
purification is needed (see Optional Step below).

### Optional Step: Density Gradient Purification (Further Purification)

If downstream analysis involves proteomics, separate and purify using a
sucrose solution or iodixanol. This step is generally not required for
functional validation assays. This procedure is documented separately as
Protocol [EV-Purification-001]{.underline} (Density Gradient
Purification).

# Critical Control Points

  ------------------------------------------------------------------------
  **⚠ Step**         **Why it is critical**    **How to avoid mistakes**
  ------------------ ------------------------- ---------------------------
  Step 4.1--4.2 ---  Tubes are thin-wall       Fill the tube close to its
  Loading and        plastic with low          capacity (e.g., top up to
  filling the        resistance to             37--38 mL with PBS for a
  ultracentrifuge    compression; if the       \~30 mL sample); verify the
  tube               liquid volume is not      volume before loading
                     nearly full, the vacuum   
                     drawn during the run will 
                     collapse/destroy the tube 
                     (see Figure 2)            

  Step 1.2 / Step    Imbalance can damage the  Use a scale or match liquid
  4.3 --- Balancing  rotor and compromise run  meniscus for centrifuge
  bottles/tubes      safety                    bottles; balance
                                               ultracentrifuge tubes
                                               precisely to 0.01 g with an
                                               analytical balance; load
                                               rotor positions
                                               symmetrically
  ------------------------------------------------------------------------

# Troubleshooting

  ---------------------------------------------------------------------------
  **Problem**           **Possible Cause**     **Solution**
  --------------------- ---------------------- ------------------------------
  Ultracentrifuge tube  Combined sample + PBS  Before the run, top up the
  collapsed/deformed    volume in the tube was liquid volume close to the
  after the run (Step   insufficient; the      tube\'s capacity (e.g., to
  4)                    resulting air gap      \~37--38 mL for a 38.5 mL
                        allows the vacuum      tube); double-check the volume
                        drawn during the run   before loading; ensure precise
                        to crush the thin-wall balancing
                        tube                   

  No visible pellet     Normal in most runs    Do not discard the tube.
  after                 --- bacterial EV       Resuspend the bottom of the
  ultracentrifugation   pellets are small and  tube as described in Step 4.7
                        translucent and are    and confirm recovery by NTA
                        often not visible by   before concluding that the run
                        eye                    failed

  TFF flux drops        Membrane fouling by    Pause and rinse or backflush
  sharply during        residual cells and     the cassette with DI water;
  filtration or         debris                 confirm that the Step 1
  concentration                                centrifugation removed the
                                               cell pellet completely;
                                               replace the cassette if the
                                               flux does not recover

  Filtrate still turbid Membrane damaged or    Re-filter through a fresh 0.2
  after 0.2 µm          bypassed through a     µm membrane; check the seal of
  filtration            poor seal; incomplete  the filter holder; check the
                        bacterial removal in   Step 1 supernatant for
                        Step 1                 carry-over of the pellet
  ---------------------------------------------------------------------------

![](./images/media/image2.jpg){width="1.3541666666666667in"
height="1.8125in"} ![](./images/media/image3.jpg){width="1.71875in"
height="2.2916666666666665in"}
![](./images/media/image4.jpg){width="1.71875in"
height="2.2916666666666665in"}

**Figure 2. Vacuum-induced collapse/deformation of the ultracentrifuge
tube**

*Left to right: (1) insufficient liquid volume in the tube (\"the volume
is not enough\"); (2) a visibly deformed tube inside the rotor bucket;
(3) close-up of the crushed tube.*

# Expected Outcome

A successful preparation yields a translucent to pale-yellow EV pellet
resuspended in \~1 mL of pre-cooled PBS, with no visible bacterial
contamination, suitable for downstream functional validation assays.
Further purification by density gradient can be performed if proteomic
analysis is needed.

# Storage and Stability

  ----------------------------------------------------------------------------
  **Item**                    **Storage Condition**   **Stability**
  --------------------------- ----------------------- ------------------------
  Concentrated sample         4°C                     1--2 days
  (pre-ultracentrifugation)                           

  Final EV                    4°C (short term) or     4°C --- short-term
  pellet/resuspension         −20°C/−80°C (long term) storage (\~1 week).
                                                      −20°C or −80°C ---
                                                      long-term storage;
                                                      aliquot before freezing
                                                      and avoid repeated
                                                      freeze--thaw cycles.
                                                      Choose the temperature
                                                      based on the intended
                                                      downstream experiment

  TFF system (idle/storage)   Internally circulated   Inhibits bacterial
                              with 1 L of 0.01 M NaOH growth

  Concentrator filter inserts Stored in 0.01 M NaOH   ---
                              or 75% ethanol          
  ----------------------------------------------------------------------------

# References

### SOPs

- EV-Purification-001 --- Density Gradient Purification (Optional Step)

### Publications

Su, Z. Z., Hou, Y., Lin, J. J., et al. (2026). *Extracellular vesicles
from biofilm and planktonic Pseudomonas aeruginosa: proteomic profiles,
iron chelation and functional implications*. *Archives of Microbiology*,
208, 239.
[[https://doi.org/10.1007/s00203-026-04794-w]{.underline}](https://doi.org/10.1007/s00203-026-04794-w)
--- detailed methods are described in this publication.

Zhu, L.-T., Huang, H.-N., Avellán-Llaguno, R. D., Qin, Y., An, X.-L.,
Su, J.-Q., Huang, Q., & Zhu, Y.-G. (2022). *Diverse functional genes
harboured in extracellular vesicles from environmental and human
microbiota*. *Journal of Extracellular Vesicles*, 11, e12292.
[[https://doi.org/10.1002/jev2.12292]{.underline}](https://doi.org/10.1002/jev2.12292)

### Instrument manuals

- *Beckman Coulter SW 32 Ti rotor manual*

- *TFF system manual*

- *Amicon® Ultra-15 concentrator tube instructions*

# Ordering Information

Consumables used in this protocol and where they are ordered.

**Table 5. Consumables and ordering information**

+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| **Item**        | **Specification** | **Supplier / Cat. No.**                                                                                                          | **Price**  |
+=================+===================+==================================================================================================================================+============+
| Ultracentrifuge | Open-Top Thinwall | Beckman Coulter, Cat. No. [[326823]{.underline}](https://www.beckman.com/supplies/tubes-and-bottles/open-top-tubes/326823)       | 326823: on |
| tubes           | PP, 38.5 mL, 25 × | (50/pack, non-sterile)                                                                                                           | request    |
|                 | 89 mm             |                                                                                                                                  |            |
|                 |                   | Sterile version: Beckman C14301 = Fisher                                                                                         | C14301:    |
|                 |                   | [[NC2259182]{.underline}](https://www.fishersci.com/shop/products/strl-38-5ml-cert-tube-48-pk/NC2259182) (48/pack)               | US\$       |
|                 |                   |                                                                                                                                  | 567.12 /   |
|                 |                   |                                                                                                                                  | 48         |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Concentrator    | Amicon® Ultra-15  | Merck Millipore, [[UFC9100]{.underline}](https://www.sigmaaldrich.com/DE/de/product/mm/ufc9100) --- pack sizes 8 / 24 / 96       | On request |
| tubes           | centrifugal       |                                                                                                                                  | (sign-in)  |
|                 | filter, 100 kDa   | 24-pack = UFC910024                                                                                                              |            |
|                 | MWCO              |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Syringe filters | Millex®-GP, 0.22  | Merck/Sigma-Aldrich, Cat. No. [[SLGPR33RS]{.underline}](https://www.sigmaaldrich.com/DE/de/product/mm/slgpr33rs)                 | € 239.00 / |
|                 | µm PES, sterile,  |                                                                                                                                  | 50         |
|                 | 33 mm             |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Filter membrane | PES, 47 mm        | Cytiva (Pall) Supor®-800, Cat. No. 60110 --- Fisher                                                                              | On request |
| --- 0.8 µm      |                   | [[17184481]{.underline}](https://www.fishersci.de/shop/products/supor-800-membrane-disc-filters-7/17184481) (100/pack)           |            |
| (pre-filter)    |                   |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Filter membrane | PES, 47 mm        | Cytiva (Pall) Supor®-450, Cat. No. 60173 --- Fisher                                                                              | ≈ SEK      |
| --- 0.45 µm     |                   | [[17134491]{.underline}](https://www.fishersci.de/shop/products/supor-450-membrane-disc-filters-9/17134491) (100/pack)           | 3,125 /    |
|                 |                   |                                                                                                                                  | 100        |
|                 |                   | Alternative: Sartorius type 15406, 15406\--47\-\-\-\-\--N                                                                        | (Fisher SE |
|                 |                   |                                                                                                                                  | list; DE   |
|                 |                   |                                                                                                                                  | price on   |
|                 |                   |                                                                                                                                  | request)   |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Filter membrane | PES, 47 mm        | Cytiva (Pall) Supor®-200, Cat. No. 60301 --- Fisher                                                                              | On request |
| --- 0.2 µm      |                   | [[17164501]{.underline}](https://www.fishersci.de/shop/products/supor-200-membrane-disc-filters-7/17164501) (100/pack)           |            |
|                 |                   |                                                                                                                                  |            |
|                 |                   | Alternative: Sartorius type 15407, 15407\--47\-\-\--MIN                                                                          |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Syringes        | 20 mL / 50 mL,    | Any supplier                                                                                                                     | ---        |
|                 | sterile           |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| 1 L centrifuge  | for SORVALL LYNX  | Thermo Scientific                                                                                                                | ---        |
| bottles         | 6000              |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Centrifuge      | 50 mL             | Any supplier                                                                                                                     | ---        |
| tubes           |                   |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| EP tubes        | 1.5 mL / 2 mL     | Eppendorf [[Safe-Lock                                                                                                            | ---        |
|                 |                   | Tubes]{.underline}](https://www.eppendorf.com/us-en/Products/Lab-Consumables/Lab-Tubes/Eppi-Eppendorf-Safe-Lock-Tubes-p-PF-8863) |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| TFF cassette    | Hydrosart®, 0.2   | Sartorius Sartocon® Slice, Order No.                                                                                             | On request |
| (Step 2 --- 0.2 | µm, 0.1 m²        | [[3051860701W]{.underline}](https://shop.sartorius.com/de/p/sartocon-slice-cassette-hydrosart-02m-01m/3051860701W--SW)           |            |
| µm filtration)  |                   |                                                                                                                                  | Around     |
|                 |                   |                                                                                                                                  | 1,500 Euro |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| TFF cassette    | Pellicon® 2 Mini  | Merck: holder [[XX42PMINI]{.underline}](https://www.sigmaaldrich.com/DE/de/product/mm/xx42pmini), module                         | Merck: on  |
| (Step 3, Method | ultrafiltration   | [[P2C100C01]{.underline}](https://www.sigmaaldrich.com/DE/de/product/mm/p2c100c01)                                               | request    |
| B --- 100 kDa)  | module, 0.1 m²,   |                                                                                                                                  |            |
|                 | 100 kDa NMWCO,    | Sartorius alternative: Sartocon® Hydrosart® 100 kDa, 0.1 m², E-Screen,                                                           | Sartorius: |
|                 | with Pellicon®    | [[3051446801E]{.underline}](https://shop.sartorius.com/de/p/sartocon-cassette-hydrosart-100kda-01m-e-screen/3051446801E--SW)     | US\$ 1,610 |
|                 | Mini cassette     |                                                                                                                                  |            |
|                 | holder            |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Pump (Step 3,   | Pump head YZ15;   | Golander [[BT100S variable-speed peristaltic                                                                                     | ---        |
| Method B)       | takes #16 and #17 | pump]{.underline}](https://www.golander.de/products/bt100s-basic-variable-speed-peristaltic-pump), 1--4 channels, 0.0001--720    |            |
|                 | tubing            | mL/min                                                                                                                           |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Wofasteril      | peracetic acid;   | Any supplier                                                                                                                     | ---        |
|                 | used as a 2%      |                                                                                                                                  |            |
|                 | (v/v) working     |                                                                                                                                  |            |
|                 | solution          |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| PBS salts       | NaCl, KCl,        | Any supplier                                                                                                                     | ---        |
|                 | Na₂HPO₄           |                                                                                                                                  |            |
|                 | (anhydrous),      |                                                                                                                                  |            |
|                 | KH₂PO₄ ---        |                                                                                                                                  |            |
|                 | analytical grade  |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| NaOH            | pellets,          | Any supplier                                                                                                                     | ---        |
|                 | analytical grade  |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+
| Ethanol         | 70--75% (v/v),    | Any supplier                                                                                                                     | ---        |
|                 | lab-grade         |                                                                                                                                  |            |
+-----------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------+

*List prices checked on 2026-10-22, excluding VAT. Prices depend on the
framework contract and on the pack size ordered --- confirm them before
ordering.*

# Revision History

  -----------------------------------------------------------------------------
  **Version**   **Date**     **Changes**                        **Author**
  ------------- ------------ ---------------------------------- ---------------
  1.0           2026-10-21   Initial version                    Haining Huang

  -----------------------------------------------------------------------------
