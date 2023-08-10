QSCORES = document.getElementById('qscores');
data = document.querySelector("meta[name='qscores_y']").content;
console.log(data);

function createQScoresPlot(data) {
    var layout = {
      width: 500,
      height: 350,
      margin: {
          l: 95,
          r: 80,
          b: 90,
          t: 10
      },
      paper_bgcolor:'rgba(0,0,0,0)',
      plot_bgcolor:'rgba(0,0,0,0)',
      modebar: {
          // vertical modebar button layout
          orientation: 'v',
          // for demonstration purposes
          bgcolor: 'white',
          color: '#a7a8a8',
          activecolor: '#9ED3CD'
        },
      
        // move the legend inside the plot on the top right
        xaxis: {
          title: {
            text: 'Q Score',
            font: {
              family: 'Rubik',
              size: 18,
              color: '#7f7f7f'
            }
          },
        },
        yaxis: {
          title: {
            text: 'Number of Reads',
            font: {
              family: 'Rubik',
              size: 18,
              color: '#7f7f7f'
            }
          }
        }
      
    }

    var config = {
      displayModeBar: true // or 'true' for always visible
    };
 
    var data = data; // string
    data = data.replace(/["']/g, "");
    y_values = JSON.parse(data);

    Plotly.plot(QSCORES, [{
    type: 'scatter',
    x: [5, 7, 10, 12, 15],
    y: y_values
    }], layout, config)
}

createQScoresPlot(data);
QSCORES.style.zIndex="100";