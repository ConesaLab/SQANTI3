const overlay = document.querySelector(".overlay");
const btnCloseModal = document.querySelector(".close-modal");
const btnOpenModal = document.querySelector(".show-modal");

const openModal = function (e) {
  if (e && e.preventDefault) e.preventDefault();
  overlay.classList.remove("hidden");
};

const closeModal = function () {
  overlay.classList.add("hidden");
};

if (btnOpenModal) btnOpenModal.addEventListener("click", openModal);
if (btnCloseModal) btnCloseModal.addEventListener("click", closeModal);


document.addEventListener("keydown", function (e) {
  if (e.key === "Escape" && !overlay.classList.contains("hidden")) {
    closeModal();
  }
});
