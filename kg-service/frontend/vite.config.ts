import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Served under /console in production (StaticFiles mount in the FastAPI app).
// In dev, /api is proxied to the backend so the SPA and API share an origin.
//
// The proxy target and dev port are configurable via env so a non-default
// backend port "just works":
//   VITE_API_TARGET=http://localhost:8001 VITE_PORT=5174 npm run dev
// In production (backend serves the built app) the SPA uses relative /api URLs,
// so it always hits whatever port the backend runs on — no config needed.
const API_TARGET = process.env.VITE_API_TARGET || "http://127.0.0.1:8000";
const DEV_PORT = Number(process.env.VITE_PORT || 5173);

export default defineConfig({
  base: "/console/",
  plugins: [react()],
  server: {
    host: "127.0.0.1",
    port: DEV_PORT,
    proxy: {
      "/api": {
        target: API_TARGET,
        changeOrigin: true,
      },
    },
  },
});
