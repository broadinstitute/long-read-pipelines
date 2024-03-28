VERSION = 0.0.1
TAG1 = us.gcr.io/broad-dsp-lrma/lr-plot_sv_metrics:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-plot_sv_metrics:latest

all: build

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
