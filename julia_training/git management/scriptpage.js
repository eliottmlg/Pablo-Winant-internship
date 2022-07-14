// Select the button
const btn = document.querySelector(".btn-toggle");
// Select the theme preference from localStorage
const currentTheme = localStorage.getItem("theme");

// If the current theme in localStorage is "dark"...
if (currentTheme == "dark") {
  // ...then use the .dark-theme class
  document.body.classList.add("dark");
}

// Listen for a click on the button 
btn.addEventListener("click", function() {
  // Toggle the .dark-theme class on each click
  document.body.classList.toggle("dark");
  
  // Let's say the theme is equal to light
  let theme = "light";
  // If the body contains the .dark-theme class...
  if (document.body.classList.contains("dark")) {
    // ...then let's make the theme dark
    theme = "dark";
  }
  // Then save the choice in localStorage
  localStorage.setItem("theme", theme);
});

var button = document.getElementById("audio-button");
var audio = document.getElementById("player");

button.addEventListener("click", function(){
  if(audio.paused){
    audio.play();
    button.innerHTML = "OFF";
  } else {
    audio.pause();
    button.innerHTML = "ON";
  }
});