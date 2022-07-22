VERSION = throwaway
TAG0 = us.gcr.io/broad-dsp-lrma/lr-gatk:$(VERSION)

all: build push

build:
	docker build -t $(TAG0) .

push:
	docker push $(TAG0)
