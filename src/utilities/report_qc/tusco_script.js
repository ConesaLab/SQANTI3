document.addEventListener("DOMContentLoaded", function () {
  // Open modal event
  document.querySelectorAll(".open-modal").forEach(function (button) {
    button.addEventListener("click", function () {
      var modalId = this.getAttribute("data-modal");
      document.getElementById(modalId).style.display = "block";
    });
  });

  // Close modal event
  document.querySelectorAll(".close-modal").forEach(function (span) {
    span.addEventListener("click", function () {
      var modalId = this.getAttribute("data-modal");
      document.getElementById(modalId).style.display = "none";
    });
  });

  // Close modal when clicking outside the modal content
  window.addEventListener("click", function (event) {
    document.querySelectorAll(".modal").forEach(function (modal) {
      if (event.target == modal) {
        modal.style.display = "none";
      }
    });
  });
});
