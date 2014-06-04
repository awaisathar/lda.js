lda.js
======

[LDA-Based Topic Modelling in Javascript](http://awaisathar.github.io/lda.js/)

Topic modelling means detecting “abstract” topics from a collection of text documents. The most common text book technique to do that is using Latent Dirichlet Allocation. Simply put, LDA is a statistical algorithm which takes documents as input and produces a list of topics. One catch is that you have to tell it how many topics you want. There’s much more to it but since this is not a tutorial post, I will stop here. (If you are interested in how it works, read the references given on the wiki page.)

![Output PNG](http://chaoticity.com/images/lda.js.png) 

Here's a Javascript version of LDA, based on my no-longer-functioning [earlier work](http://chaoticity.com/lda-based-topic-modelling-in-javascript/). For testing, I use a subset of the SMS Spam Corpus [available here](http://www.dt.fee.unicamp.br/~tiago/smsspamcollection/) (and thus take no responsibility of the inappropriateness of the text within :) ). Each topic is represented as a word cloud; the larger a word, the more weight it has in the topic. The source sentences are displayed again with a bar which shows the percentage distribution of topics for that sentence. Hovering on each area in the bar would show you the words in the topic. You can of course replace it with any other text, change the number of topics using the slider, and press the 'Analyse' button to see it work.
