function bindTuscoModalHandlers() {
  // Open modal event
  document.querySelectorAll('.open-modal').forEach(function (button) {
    if (!button.__tuscoBound) {
      button.addEventListener('click', function () {
        var modalId = this.getAttribute('data-modal');
        var el = document.getElementById(modalId);
        if (el) el.style.display = 'block';
      });
      button.__tuscoBound = true;
    }
  });

  // Close modal event
  document.querySelectorAll('.close-modal').forEach(function (span) {
    if (!span.__tuscoBound) {
      span.addEventListener('click', function () {
        var modalId = this.getAttribute('data-modal');
        var el = document.getElementById(modalId);
        if (el) el.style.display = 'none';
      });
      span.__tuscoBound = true;
    }
  });

  // Close modal when clicking outside the modal content
  if (!window.__tuscoWindowBound) {
    window.addEventListener('click', function (event) {
      document.querySelectorAll('.modal').forEach(function (modal) {
        if (event.target === modal) {
          modal.style.display = 'none';
        }
      });
    });
    window.__tuscoWindowBound = true;
  }
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', bindTuscoModalHandlers);
} else {
  bindTuscoModalHandlers();
}
