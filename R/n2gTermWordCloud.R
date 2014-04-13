# This function uses output from n2gExploreTerm.py (the raw text) to perform
# stop word removal, stemming, and generation of a word cloud.  The goal of this
# is to understand the "gist" of the papers that the terms are derived from.

# Load libraries
library('tm')
library('wordcloud')

# Read in raw data files - make sure your raw text file is in this folder
rawtext = '/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/nlp'

# Create a document "corpus" (we only have one document)
corpy = Corpus(DirSource(rawtext),readerControl = list(language="en"))

# Do inspect(corpy) to see changes
# Get rid of white space
corpy = tm_map(corpy, stripWhitespace)
# We already did .lower() in python, but just in case
corpy = tm_map(corpy, tolower)
# Remove English stop words
corpy = tm_map(corpy, removeWords, stopwords("english"))
# Remove punctuation
corpy = tm_map(corpy, removePunctuation)
# Word stemming
corpy = tm_map(corpy, stemDocument)

# MAKE WORDCLOUD!
wordcloud(corpy, scale=c(5,0.5), max.words=100, random.order=FALSE, rot.per=0.35, use.r.layout=FALSE, colors=brewer.pal(8,"Dark2"))