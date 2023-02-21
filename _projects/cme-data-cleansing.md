---
layout: page
title: CME Group
description: Cleansing of financial market data - a generic framework
img: assets/img/projects/cme-data-cleansing/cover.png
importance: 2
category: work
---

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/projects/cme-data-cleansing/framework.png" title="State-cluster-behaviour modelling" class="img-fluid rounded z-depth-1" %}
    </div>
</div>

---

CME Clearing (referred to as "CME") uses a historical simulation approach to generate risk scenarios for client's portfolio. The distribution of profits and losses is constructed by taking the current portfolio, and subjecting it to the actual changes in the risk factors experienced during each day of a historical period. It is crucial to ensure the historical data to be cleansed before entering into the simulation engine.

Our project aims to construct a generic framework that is capable of consistently describing the mainstream methods and processes that are currently in use at CME global quant team to cleanse financial data associated with a variety of products and asset classes. The framework provides a quick guidance on the workow to follow and an encyclopedia of techniques to use for a data cleansing task. This will benefit CME by (1) accelerating risk modelling research and promoting fast launch of new products, and (2) that by comparing different methods and processes in the framework, CME could seek for methodology alignment across the global team.

We start from defining the concept _financial pricing unit_ where the data cleansing task is carried out on. A financial pricing unit (FPU) is a structure of a financial variable used as input by a pricer to evaluate any financial product. Next we summarise all FPUs needed for pricing asset classes and products in scope. These products include linear and non-linear derivatives on foreign exchange rate, interest rate and commodity (energy). By arguing for the nite-dimensionality and discreteness properties of the FPU data, we come up with a mathematical representation for the data to be studied.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/projects/cme-data-cleansing/fpu.png" title="Pipeline" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Map of financial pricing units required for pricing different products in scope.
</div>


With data representation being dened, we move to propose the principles that clarify what properties the cleansed data should have, such as completeness, outlier-free, no violation to pricing constraints, and so on. The principles give rise to three key processes to cleanse data, namely price validation, outlier detection and data completion. These processes are organised in a particular order to form the generic framework. There are possibly multiple approaches to accomplish a process. We call each approach a method. We then explain various methods in details.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/projects/cme-data-cleansing/process.png" title="Pipeline" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    A generic framework of data cleansing processes.
</div>

Readers who are interested in this project could refer to the [lay report](https://www.maths.ox.ac.uk/system/files/attachments/CME-Victor-LayReport_Final.pdf) and are very welcome to contact me through email for more details.
